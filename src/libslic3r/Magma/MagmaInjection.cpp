#include "MagmaInjection.hpp"
#include "MagmaTubeMap.hpp"

#include "../GCode.hpp"
#include "../GCode/GCodeProcessor.hpp"
#include "../MultiPoint.hpp"
#include "../ShortestPath.hpp"
#include "../I18N.hpp"
#include "../format.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

#include <boost/log/trivial.hpp>

namespace Slic3r {
namespace magma {

// One renderable column of an injection manifold: a centerline with a per-point render
// width (parallel vectors), so the column can taper layer-by-layer where the part clips it.
struct TubeVizLine { std::vector<Vec3d> pts; std::vector<double> widths; };

// Keep at most max_n points per column (endpoints + evenly spaced), preserving widths. Not RDP:
// the preview animates the fill point by point, so a straight column still needs midpoints.
static void subsample_with_widths(std::vector<Vec3d>& pts, std::vector<double>& widths, size_t max_n)
{
    if (pts.size() <= max_n || max_n < 2)
        return;
    std::vector<Vec3d>  np;  np.reserve(max_n);
    std::vector<double> nw;  nw.reserve(max_n);
    for (size_t i = 0; i < max_n; ++i) {
        const size_t idx = (i * (pts.size() - 1)) / (max_n - 1); // includes first and last
        np.push_back(pts[idx]);
        nw.push_back(idx < widths.size() ? widths[idx] : 0.4);
    }
    pts = std::move(np); widths = std::move(nw);
}

// Format the manifold as one MAGMA_TUBE comment. lines[0] is the hub; the rest branch from
// its last point. Points are x,y,z,width (parsed in GCodeProcessor::process_tags).
static std::string format_tube_viz_comment(const std::vector<TubeVizLine>& lines)
{
    std::ostringstream oss;
    oss << ";" << GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Magma_Tube)
        << "np=" << lines.size() << " n=";
    for (size_t i = 0; i < lines.size(); ++i) { if (i) oss << ','; oss << lines[i].pts.size(); }
    oss << " pts=";
    bool first = true;
    for (const TubeVizLine& line : lines)
        for (size_t k = 0; k < line.pts.size(); ++k) {
            if (!first) oss << ';';
            first = false;
            char buf[80];
            snprintf(buf, sizeof(buf), "%.3f,%.3f,%.3f,%.3f",
                     line.pts[k].x(), line.pts[k].y(), line.pts[k].z(),
                     k < line.widths.size() ? line.widths[k] : 0.1);
            oss << buf;
        }
    oss << '\n';
    return oss.str();
}

std::string park_and_set_temp(
    GCode& gcodegen,
    bool park_enabled,
    double print_z,
    double park_z_hop,
    double extra_retract,
    int target_temp,
    const char* z_comment,
    const char* xy_comment)
{
    if (park_enabled) {
        ParkResult park = gcodegen.safe_park().find_safe_position(
            gcodegen.layer(), nullptr, gcodegen.last_pos());
        return SafeParkPosition::park_and_set_temp(
            gcodegen, park, print_z, park_z_hop, extra_retract,
            target_temp, z_comment, xy_comment);
    }
    // Not parking — just retract and set temperature.
    std::string gcode;
    gcode += gcodegen.retract(false, false);
    gcode += gcodegen.writer().set_temperature(target_temp, true);
    return gcode;
}

std::vector<std::string> InjectionDiagnostics::messages() const
{
    std::vector<std::string> out;
    if (tubes_no_volume > 0)
        out.push_back(Slic3r::format(
            L("%1% Magma U-tube(s) were not injected: no injectable cavity was measured in the "
              "printed lattice, so those channels are printing hollow and unreinforced. This is "
              "usually a tube clipped away to nothing at a part boundary, or a tube height that "
              "leaves no cavity between the window and the cap."),
            tubes_no_volume));
    if (layers_no_fill_factor > 0)
        out.push_back(Slic3r::format(
            L("Every Magma injection on %1% layer(s) was dropped because Tube fill factor is not "
              "greater than zero. The lattice on those layers prints hollow."),
            layers_no_fill_factor));
    if (injections_speed_floored > 0)
        out.push_back(Slic3r::format(
            L("%1% Magma injection(s) ran at 1 mm\u00b3/s instead of the slower rate configured: "
              "injection speed is floored at 1 mm\u00b3/s. Rates below that hold the nozzle in "
              "the cell long enough to melt it. Raise Magma injection speed to at least "
              "1 mm\u00b3/s to choose the rate yourself."),
            injections_speed_floored));
    if (injections_invented_speed > 0)
        out.push_back(Slic3r::format(
            L("%1% Magma injection(s) ran at 1 mm\u00b3/s, a rate nobody chose: injection speed is "
              "set to 0 (\"use the filament's max volumetric speed\") but the filament declares no "
              "max volumetric speed. Set one or the other."),
            injections_invented_speed));
    if (injections_unlabelled > 0)
        out.push_back(Slic3r::format(
            L("%1% Magma injection(s) are not inside an object-exclusion block, so they will still "
              "run if you cancel that object mid-print \u2014 the nozzle will press into and extrude "
              "onto a part that is no longer being printed."),
            injections_unlabelled));
    return out;
}

std::string generate_injection_gcode(
    GCode& gcodegen,
    const MagmaTubeMap& tube_map,
    const std::vector<InjectionPoint>& points,
    double print_z,
    InjectionDiagnostics& diag)
{
    if (points.empty())
        return {};

    std::string gcode;
    const auto& config = gcodegen.config();
    char buf[256];

    // Every Z below is a nozzle Z. The layer prints at print_z + z_offset (GCode::change_layer),
    // so the seal, plunge and iron are measured from there, and the MAGMA_TUBE comment carries
    // the same absolute Z that GCodeProcessor subtracts z_offset from.
    const double z_offset = config.z_offset.value;
    const double layer_z  = print_z + z_offset;

    // --- Config values ---
    double injection_speed_vol = config.magma_injection_speed.value;
    double fill_factor       = config.magma_tube_fill_factor.value;
    if (fill_factor <= 0) {
        // Unreachable from the GUI (min 0.1); possible from a hand-edited or foreign 3MF.
        BOOST_LOG_TRIVIAL(error) << "Magma: tube fill factor is " << fill_factor
                                 << " (must be > 0); skipping all injections on this layer.";
        ++diag.layers_no_fill_factor;
        return {};
    }

    int    dwell_ms          = config.magma_injection_dwell.value;
    bool inj_retract         = config.magma_injection_retract.value;

    // Crater wipe ("plow the displaced rim back, scrape nozzle clean").
    bool   wipe_enabled      = config.magma_injection_iron.value;
    int    wipe_turns_cfg    = config.magma_injection_iron_turns.value;  // 0 = auto
    const double wipe_speed_mms = resolve_iron_speed(config.magma_injection_iron_speed.value,
                                                     config.ironing_speed.value);

    // Get extruder info (ToolOrdering handles filament switching before we get here)
    unsigned int extruder_id = gcodegen.writer().filament()->id();
    double filament_diameter = config.filament_diameter.get_at(extruder_id);
    double filament_area = (PI / 4.0) * filament_diameter * filament_diameter;

    // Injection tip. The seal and plunge depths are solved per tube below; the maths lives in
    // MagmaTriangleCell.hpp, shared with Print::validate() and the injection timing model.
    const double seal_flat  = resolve_nozzle_flat(config.magma_nozzle_outer_diameter.get_at(extruder_id),
                                                  config.nozzle_diameter.get_at(extruder_id));
    const double cone_deg   = config.magma_nozzle_cone_half_angle.get_at(extruder_id);
    const double seal_press = std::max(0.0, config.magma_seal_press.value);
    // Plunge: descend a further plunge_depth below seal_depth during the injection, so the hot
    // tip sinks into the softening tube top and keeps the opening sealed as the channel fills.
    const double plunge_cfg = requested_plunge_depth(config.magma_injection_plunge.value,
                                                     config.magma_injection_plunge_depth.value);

    // Raw volume→E conversion: 1/cross_section, without filament_flow_ratio.
    // Injection volume is geometrically computed from tube dimensions; applying
    // flow ratio would double-correct (fill_factor already controls fill amount).
    double e_per_mm3 = 1.0 / filament_area;

    // Injection volumetric speed. Filling at the material's melt limit (0, the default) fills
    // more reliably than a fixed rate: the tube is narrow and the melt cools the whole way down.
    // Print::validate rejects 0 on a filament with no max volumetric speed, so `invented` means
    // the config bypassed validation (hand-edited or foreign 3MF); both substitutions are
    // reported through diag.
    const InjectionRate rate = resolve_injection_rate(
        injection_speed_vol, config.filament_max_volumetric_speed.get_at(extruder_id));
    if (rate.invented) {
        BOOST_LOG_TRIVIAL(error)
            << "Magma: injection speed is 0 (use filament max volumetric speed) but filament "
            << extruder_id << " declares no max volumetric speed; injecting at 1 mm3/s. "
            << "Set either Magma injection speed or the filament's max volumetric speed.";
        diag.injections_invented_speed += (int) points.size();
    }
    if (rate.floored)
        diag.injections_speed_floored += (int) points.size();
    const double vol_speed = rate.mm3_s;

    // Convert volumetric speed to filament feedrate
    double feedrate_mms = vol_speed / filament_area;
    double feedrate_mmmin = feedrate_mms * 60.0;

    // Z speed for the seal and lift moves: the injecting extruder's Z travel speed, with the
    // same fallback as GCodeWriter::_travel_to_z. A slow Z lingers on the hot tube top.
    double z_speed_mms = config.travel_speed_z.get_at(extruder_id);
    if (z_speed_mms <= 0.)
        z_speed_mms = config.travel_speed.get_at(extruder_id);
    int z_feedrate = (int)(z_speed_mms * 60.0);
    if (z_feedrate < 60) z_feedrate = 60;

    // --- Injection loop ---
    float display_dim = tube_map.interior_width();

    for (const auto& pt : points) {
        double volume = pt.volume_mm3 * fill_factor;
        if (volume <= 0) {
            BOOST_LOG_TRIVIAL(warning)
                << "Magma: tube pair " << pt.pair_index << " has no measured cavity volume ("
                << pt.volume_mm3 << " mm3); skipping its injection.";
            ++diag.tubes_no_volume;
            continue;
        }

        // Seal to this tube's measured cap opening: boundary-clipped caps are smaller and the
        // nominal depth would over-press them. The plunge is re-clamped to the new seal depth.
        const double seal_depth   = auto_seal_depth(tube_map.cap_opening_diameter(tube_map.u_tube_pairs()[pt.pair_index]),
                                                    seal_flat, cone_deg, seal_press);
        const double plunge_depth = clamp_plunge_depth(seal_depth, plunge_cfg);

        // Travel to injection point. The built-in travel_to() handles retraction,
        // avoid-crossing-perimeters, and the printer's own z-hop; the unretract
        // below also undoes any lift it applied.
        Point scaled_pos(scale_(pt.position.x()), scale_(pt.position.y()));
        gcode += gcodegen.travel_to(scaled_pos, erMagmaInjection, "move to injection point");

        // Unretract — undoes the previous injection's retract (state-tracked via
        // GCodeWriter::retract) and any travel retract/lift from travel_to().
        gcode += gcodegen.unretract();

        // Seal: one fast move down into the surface, sealing against the tube opening
        if (seal_depth > 0) {
            sprintf(buf, "G1 Z%.3f F%d ; magma seal\n", layer_z - seal_depth, z_feedrate);
            gcode += buf;
        }

        // Set role and display dimensions for this injection.
        // Emitted per-injection (after unretract) so that only the actual
        // injection extrusion is classified as erMagmaInjection.
        sprintf(buf, ";%s%s\n",
                GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Role).c_str(),
                ExtrusionEntity::role_to_string(erMagmaInjection).c_str());
        gcode += buf;
        sprintf(buf, ";%s%g\n",
                GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Height).c_str(),
                display_dim);
        gcode += buf;
        sprintf(buf, ";%s%g\n",
                GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Width).c_str(),
                display_dim);
        gcode += buf;

        // Preview geometry: one MAGMA_TUBE comment carrying the manifold -- a hub column (cap
        // down to the window) plus one leg per vent (cell_b + extra_vents) branching from the
        // hub's bottom. Columns follow the per-layer cell centres and the per-layer bore.
        const auto&  pair    = tube_map.u_tube_pairs()[pt.pair_index];
        const int    cap     = pair.pair_end_layer;
        const int    wc      = pt.window_center_layer;
        const double iw      = tube_map.interior_width();
        const double z_floor = tube_map.print_z(pair.pair_start_layer) + iw * 0.5;
        auto z_at = [&](int L) { return std::max(tube_map.print_z(L), z_floor) + z_offset; };

        std::vector<TubeVizLine> lines;
        TubeVizLine hub;
        for (int L = cap; L >= wc; --L) {
            // The cell centre moves per layer with the spiral interlock; the cap point is the
            // injection XY, where the nozzle sits.
            const Vec2d hub_xy = (L == cap) ? pair.injection_center
                                            : tube_map.lattice_at(L).cell_center(pair.cell_a);
            hub.pts.push_back({ hub_xy.x(), hub_xy.y(),
                                (L == cap) ? layer_z : z_at(L) });
            hub.widths.push_back(tube_map.cell_bore_at(pair.cell_a, L));
        }
        subsample_with_widths(hub.pts, hub.widths, 16);
        const Vec3d junction = hub.pts.empty()
            ? Vec3d(pair.injection_center.x(), pair.injection_center.y(), z_at(wc))
            : hub.pts.back();
        lines.push_back(std::move(hub));
        // Each leg runs from the junction across the window to its vent, then up the vent.
        auto add_leg = [&](const magma::CellId& vent) {
            TubeVizLine leg;
            // The junction point takes the hub bore, so the connector tapers from the hub
            // mouth down to the vent bore.
            leg.pts.push_back(junction);
            leg.widths.push_back(tube_map.cell_bore_at(pair.cell_a, wc));
            for (int L = wc; L <= cap; ++L) {
                Vec2d c = tube_map.lattice_at(L).cell_center(vent);
                leg.pts.push_back({ c.x(), c.y(), (L == cap) ? layer_z : z_at(L) });
                leg.widths.push_back(tube_map.cell_bore_at(vent, L));
            }
            subsample_with_widths(leg.pts, leg.widths, 16);
            if (leg.pts.size() >= 2)
                lines.push_back(std::move(leg));
        };
        add_leg(pair.cell_b);
        for (const auto& ev : pair.extra_vents)
            add_leg(ev);
        // Shift from object coordinates to the writer's: the writer sits at this injection
        // point after the travel above, and hub[0] is the object-space injection centre. The
        // preview then applies only its standard offsets, as for any move.
        const Vec2d obj_off(gcodegen.writer().get_position().x() - pair.injection_center.x(),
                            gcodegen.writer().get_position().y() - pair.injection_center.y());
        for (auto& line : lines)
            for (auto& p : line.pts) { p.x() += obj_off.x(); p.y() += obj_off.y(); }
        if (!lines[0].pts.empty())
            gcode += format_tube_viz_comment(lines);

        double filament_length = volume * e_per_mm3;
        Vec2d xy(gcodegen.writer().get_position().x(),
                 gcodegen.writer().get_position().y());

        // Inject while plunging: fold the Z descent into the extrude moves so the nozzle sinks
        // as the channel fills. F is the cartesian speed, and the cartesian distance is only
        // plunge_depth, so F is scaled for the move to last as long as the extrusion:
        //   F = vol_feedrate * (plunge_depth / filament_length)
        // An unscaled F would dump the E at the extruder's max speed and ooze around the nozzle.
        // This E:cartesian ratio trips Klipper's max_extrude_cross_section, which must be raised
        // in printer.cfg. Without plunge, injection is pure-E moves paced by the extruder.
        // Moves are split only to keep each one's E under the firmware's extrude-only limit.
        const int N = std::max(1, std::min(8, (int) std::ceil(filament_length / 40.0)));
        if (plunge_depth > 0.0) {
            const double f_inject = feedrate_mmmin * plunge_depth / std::max(1e-4, filament_length);
            gcode += gcodegen.writer().set_speed(f_inject);
            for (int k = 0; k < N; ++k) {
                double z = layer_z - (seal_depth + plunge_depth * double(k + 1) / double(N));
                gcode += gcodegen.writer().extrude_to_xyz(
                    Vec3d(xy.x(), xy.y(), z), filament_length / double(N),
                    N > 1 ? "injection segment" : "Magma injection");
            }
        } else {
            gcode += gcodegen.writer().set_speed(feedrate_mmmin);
            for (int k = 0; k < N; ++k)
                gcode += gcodegen.writer().extrude_to_xy(
                    xy, filament_length / double(N), N > 1 ? "injection segment" : "Magma injection");
        }

        // Dwell: hold nozzle sealed while plastic spreads through tube
        if (dwell_ms > 0) {
            sprintf(buf, "G4 P%d ; injection dwell\n", dwell_ms);
            gcode += buf;
        }

        // --- Finish: break the seal, retract, then wipe the crater ---
        // Lift to crack the seal before retracting, so the retract doesn't pull the fresh plug
        // back up through the sealed interface. 0.3mm relieves the contact pressure at any depth.
        const double break_lift = 0.3;
        double deep_z = layer_z - seal_depth - plunge_depth;   // nozzle depth after plunge
        sprintf(buf, "G1 Z%.3f F%d ; injection break-lift\n", deep_z + break_lift, z_feedrate);
        gcode += buf;
        if (inj_retract)
            gcode += gcodegen.writer().retract();

        if (wipe_enabled) {
            // Crater ironing: spiral the nozzle inward over the injection point at layer
            // height, so the angled cone plows the displaced rim back into the crater and
            // irons the surface flat, while scraping the nozzle clean so it doesn't string to
            // the next tube.

            // Preview classification: role Ironing, and WIPE markers so these non-extruding
            // moves show as the Wipe move type rather than under the injection role.
            sprintf(buf, ";%s%s\n",
                    GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Role).c_str(),
                    ExtrusionEntity::role_to_string(erIroning).c_str());
            gcode += buf;
            gcode += ";" + GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Wipe_Start) + "\n";

            // Centre = current nozzle XY (the injection point). All spiral points
            // are computed relative to it, in G-code (origin-applied) coordinates.
            Vec2d c(gcodegen.writer().get_position().x(),
                    gcodegen.writer().get_position().y());
            const double r_flat = seal_flat / 2.0;

            // The whole pass presses at layer height, neighbouring cells included. A neighbour's
            // open air hole is safe: at layer height the flat has no material to push into it,
            // and the spiral runs outer -> inner, so everything it displaces travels toward our
            // own crater. That only holds if the nozzle is already at layer height when it
            // leaves the centre -- see the lift below. Accepted residual risks: the hot flat can
            // slump a thin hole rim it lingers over, and the spiral's tangential component
            // smears slightly sideways.

            const CraterIron iron = plan_crater_iron(
                seal_flat, cone_deg, seal_depth + plunge_depth, config.magma_injection_iron_margin.value, wipe_turns_cfg);
            const double start_R    = iron.start_R;
            const int    wipe_turns = iron.turns;

            int    wf       = std::max(60, (int)(wipe_speed_mms * 60.0));  // mm/s -> mm/min
            const int seg_per_rev = 16;                              // circle smoothness

            // The raw G1s below bypass GCodeWriter, which is what subtracts the plate offset
            // (GCodeWriter::travel_to_xy), so subtract it here or the pass lands a plate width
            // away on any plate but the first. Z has no offset; set_position() takes pre-offset XY.
            const Vec2f plate_off = gcodegen.writer().get_xy_offset();

            // Rise straight up to layer height before moving outward. After the break-lift the
            // nozzle can still be below the surface (seal + plunge deeper than the lift), and a
            // diagonal first move would drag the flat outward through the crater wall -- the one
            // move in this pass that carries material away from our crater.
            if (deep_z + break_lift < layer_z) {
                sprintf(buf, "G1 Z%.3f F%d ; crater iron lift to surface\n", layer_z, z_feedrate);
                gcode += buf;
            }

            // Iron path: one full rim rotation to shear the whole rim evenly, then a spiral
            // inward to the centre. The radius shrinks to 0, so a fixed seg/rev over-resolves
            // the inner turns; Douglas-Peucker at 0.08mm (under the rim circle's chord error)
            // collapses them -- far fewer G-code moves and preview steps for the same path.
            Points iron_pts;
            iron_pts.reserve(static_cast<size_t>(seg_per_rev) * (wipe_turns + 1));
            auto add_iron = [&](double rad, double ang) {
                iron_pts.emplace_back(scale_(c.x() + rad * std::cos(ang)), scale_(c.y() + rad * std::sin(ang)));
            };
            for (int s = 1; s <= seg_per_rev; ++s)                          // rim
                add_iron(start_R, double(s) / double(seg_per_rev) * 2.0 * PI);
            const int total = std::max(1, wipe_turns * seg_per_rev);        // spiral inward
            for (int s = 1; s <= total; ++s) {
                const double frac = double(s) / double(total);              // 0..1, outer -> centre
                add_iron(start_R * (1.0 - frac), frac * wipe_turns * 2.0 * PI);
            }
            for (const Point& p : MultiPoint::_douglas_peucker(iron_pts, scale_(0.08))) {
                sprintf(buf, "G1 X%.3f Y%.3f Z%.3f F%d ; crater iron\n",
                        unscale<double>(p.x()) - plate_off.x(), unscale<double>(p.y()) - plate_off.y(), layer_z, wf);
                gcode += buf;
            }
            // One flat-width stroke across the centre to flatten the gathered mound.
            sprintf(buf, "G1 X%.3f Y%.3f Z%.3f F%d ; crater iron flatten\n",
                    c.x() - r_flat - plate_off.x(), c.y() - plate_off.y(), layer_z, wf); gcode += buf;
            sprintf(buf, "G1 X%.3f Y%.3f Z%.3f F%d ; crater iron flatten\n",
                    c.x() + r_flat - plate_off.x(), c.y() - plate_off.y(), layer_z, wf); gcode += buf;
            // Return to centre at layer height and resync the writer position (the
            // raw moves above bypassed its tracking) so the next travel_to plans
            // its path/avoid-crossing from the right spot.
            sprintf(buf, "G1 X%.3f Y%.3f Z%.3f F%d ; crater iron end\n",
                     c.x() - plate_off.x(), c.y() - plate_off.y(), layer_z, wf);
            gcode += buf;
            gcode += ";" + GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Wipe_End) + "\n";
            gcodegen.writer().set_position(Vec3d(c.x(), c.y(), layer_z));
        } else {
            // No wipe: just return to layer height and resync the writer's Z.
            sprintf(buf, "G1 Z%.3f F%d ; magma seal release\n", layer_z, z_feedrate);
            gcode += buf;
            Vec3d p = gcodegen.writer().get_position(); p.z() = layer_z;
            gcodegen.writer().set_position(p);
        }
        // No manual z-hop: the next iteration's travel_to() handles lift + travel.
    }

    // Reset forced dimensions
    sprintf(buf, ";%s0\n",
            GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Height).c_str());
    gcode += buf;
    sprintf(buf, ";%s0\n",
            GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Width).c_str());
    gcode += buf;

    return gcode;
}

} // namespace magma
} // namespace Slic3r
