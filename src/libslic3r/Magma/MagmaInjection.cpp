#include "MagmaInjection.hpp"
#include "MagmaTubeMap.hpp"
#include "MagmaSpiralOffset.hpp"

#include "../GCode.hpp"
#include "../GCode/GCodeProcessor.hpp"
#include "../Print.hpp"
#include "../ShortestPath.hpp"

#include <igl/ramer_douglas_peucker.h>

#include <algorithm>
#include <cmath>
#include <sstream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace Slic3r {
namespace magma {

std::vector<InjectionPoint> collect_injection_points(
    const MagmaTubeMap& tube_map,
    int layer_id)
{
    std::vector<InjectionPoint> points;

    const auto& pairs = tube_map.u_tube_pairs();
    for (int i = 0; i < static_cast<int>(pairs.size()); ++i) {
        const auto& pair = pairs[i];
        if (pair.pair_end_layer != layer_id)
            continue;
        if (pair.volume_mm3 <= 0)
            continue;

        InjectionPoint pt;
        // Use spiral-offset lattice for this layer so injection aligns with drawn pattern
        TriangleLattice lattice = lattice_for_layer(
            tube_map.cell_spacing(), tube_map.spiral_params(), layer_id);
        pt.position   = lattice.cell_center(pair.cell_a);
        pt.volume_mm3 = pair.volume_mm3;
        pt.pair_index = i;
        pt.start_layer = pair.pair_start_layer;
        // Window center: halfway through the window gap region at the tube bottom
        int wh = tube_map.window_height_layers();
        pt.window_center_layer = pair.pair_start_layer + wh / 2;
        points.push_back(pt);
    }

    // Order injection points for minimal travel using OrcaSlicer's greedy
    // TSP solver (KD-tree backed).  Much better than a raster sweep for
    // irregular point distributions near model boundaries.
    if (points.size() > 1) {
        Points scaled_pts;
        scaled_pts.reserve(points.size());
        for (const auto& pt : points)
            scaled_pts.push_back(Point(scale_(pt.position.x()), scale_(pt.position.y())));
        std::vector<size_t> order = chain_points(scaled_pts);
        std::vector<InjectionPoint> ordered;
        ordered.reserve(points.size());
        for (size_t idx : order)
            ordered.push_back(std::move(points[idx]));
        points = std::move(ordered);
    }

    return points;
}

// Generate the spiral-following center-line of a U-tube and simplify with
// Douglas-Peucker.  Returns waypoints tracing: top of cell_a → descent →
// window crossing → ascent → top of cell_b.
static std::vector<Vec3d> build_tube_viz_waypoints(
    const MagmaTubeMap& tube_map,
    const UTubePair& pair,
    double layer_z,         // Z at pair_end_layer
    int    window_center_layer)
{
    const auto& sp = tube_map.spiral_params();
    double cs = tube_map.cell_spacing();
    float  iw = tube_map.interior_width();

    // Helper: Z height for a given layer index (uses actual per-layer print_z)
    auto z_for_layer = [&](int L) -> double {
        return tube_map.print_z(L);
    };

    // Z offsets: raise tube bottom so rendered cylinder doesn't poke through floor
    double z_bot_raw = z_for_layer(pair.pair_start_layer);
    double z_bot = z_bot_raw + iw / 2.0;  // sit ON the floor
    double z_window = z_for_layer(window_center_layer);
    // Clamp window Z above the raised bottom
    if (z_window < z_bot)
        z_window = z_bot;

    std::vector<Vec3d> full_path;

    // Phase 1: Descend through cell A (top → window center)
    for (int L = pair.pair_end_layer; L >= window_center_layer; --L) {
        TriangleLattice lat = lattice_for_layer(cs, sp, L);
        Vec2d c = lat.cell_center(pair.cell_a);
        double z = (L == pair.pair_end_layer) ? layer_z :
                   (L <= window_center_layer) ? z_window :
                   std::max(z_for_layer(L), z_bot);
        full_path.push_back({c.x(), c.y(), z});
    }

    // Phase 2: Cross to cell B at window center Z
    {
        TriangleLattice lat = lattice_for_layer(cs, sp, window_center_layer);
        Vec2d c = lat.cell_center(pair.cell_b);
        full_path.push_back({c.x(), c.y(), z_window});
    }

    // Phase 3: Ascend through cell B (window center → top)
    for (int L = window_center_layer + 1; L <= pair.pair_end_layer; ++L) {
        TriangleLattice lat = lattice_for_layer(cs, sp, L);
        Vec2d c = lat.cell_center(pair.cell_b);
        double z = (L == pair.pair_end_layer) ? layer_z :
                   std::max(z_for_layer(L), z_bot);
        full_path.push_back({c.x(), c.y(), z});
    }

    // Simplify with 3D Ramer-Douglas-Peucker.
    // When spirals are off, each phase is a straight line so RDP reduces
    // ~20-120 points down to ~5 (top_A, bot_A, window_B, bot_B, top_B).
    // When spirals are on, RDP keeps points where the helix bends noticeably.
    if (full_path.size() > 2) {
        Eigen::MatrixXd P(full_path.size(), 3);
        for (size_t i = 0; i < full_path.size(); ++i)
            P.row(i) = full_path[i].transpose();

        Eigen::MatrixXd S;
        Eigen::VectorXi J;
        igl::ramer_douglas_peucker(P, double(iw) * 0.1, S, J);

        std::vector<Vec3d> simplified;
        simplified.reserve(S.rows());
        for (int i = 0; i < S.rows(); ++i)
            simplified.push_back(S.row(i).transpose());
        return simplified;
    }
    return full_path;
}

// Format waypoints as a MAGMA_TUBE G-code comment.
static std::string format_tube_viz_comment(const std::vector<Vec3d>& waypoints, float width)
{
    std::ostringstream oss;
    oss << "; MAGMA_TUBE n=" << waypoints.size() << " w=" << width << " pts=";
    for (size_t i = 0; i < waypoints.size(); ++i) {
        if (i > 0) oss << ';';
        char buf[64];
        snprintf(buf, sizeof(buf), "%.3f,%.3f,%.3f",
                 waypoints[i].x(), waypoints[i].y(), waypoints[i].z());
        oss << buf;
    }
    oss << '\n';
    return oss.str();
}

std::string generate_injection_gcode(
    GCode& gcodegen,
    const MagmaTubeMap& tube_map,
    const std::vector<InjectionPoint>& points,
    double layer_z,
    double actual_layer_height)
{
    if (points.empty())
        return {};

    std::string gcode;
    const auto& config = gcodegen.config();
    char buf[256];

    // --- Config values ---
    int injection_temp       = config.magma_injection_temp.value;
    double injection_speed_vol = config.magma_injection_speed.value;
    bool iron_tube_ends      = config.magma_iron_tube_ends.value;
    bool park_enabled        = config.magma_injection_park.value;
    int dwell_ms             = config.magma_injection_dwell.value;
    double fill_factor       = config.magma_tube_fill_factor.value;
    bool z_slam              = config.magma_injection_z_slam.value;
    int inj_filament         = config.magma_injection_filament.value;

    // --- Filament switch (before any injection) ---
    unsigned int original_filament = gcodegen.writer().filament()->id();
    bool need_filament_switch = (inj_filament > 0 &&
                                 (unsigned int)(inj_filament - 1) != original_filament);

    if (need_filament_switch) {
        gcode += "; Magma injection: switching to injection filament\n";
        gcode += gcodegen.set_extruder(inj_filament - 1, layer_z);
    }

    // Get extruder info (may have changed after filament switch)
    unsigned int extruder_id = gcodegen.writer().filament()->id();
    double filament_diameter = config.filament_diameter.get_at(extruder_id);
    double filament_area = (M_PI / 4.0) * filament_diameter * filament_diameter;

    // Flow-adjusted E-per-mm3: includes filament_flow_ratio and filament cross-section.
    // Use this for volume→E conversion so injection respects the user's flow calibration.
    // filament_area above is still needed for feedrate calc (physical, not flow-adjusted).
    double e_per_mm3 = gcodegen.writer().filament()->e_per_mm3();

    // Injection volumetric speed.
    // Default (0) uses the infill volumetric rate — this represents the
    // hotend's proven melt capacity during normal printing.
    double vol_speed;
    if (injection_speed_vol > 0) {
        vol_speed = injection_speed_vol;
    } else {
        double max_vol = config.filament_max_volumetric_speed.get_at(extruder_id);
        if (max_vol > 0) {
            vol_speed = max_vol;
        } else {
            float nozzle_d = config.nozzle_diameter.get_at(extruder_id);
            float line_width = config.sparse_infill_line_width.get_abs_value(nozzle_d);
            if (line_width <= 0) line_width = nozzle_d;
            float lh = static_cast<float>(actual_layer_height);
            double mm3_per_mm = lh * (line_width - lh * (1.0 - M_PI / 4.0));
            vol_speed = config.sparse_infill_speed.value * mm3_per_mm;
        }
    }

    // Cap at hotend melt rate
    double max_vol = config.filament_max_volumetric_speed.get_at(extruder_id);
    if (max_vol > 0)
        vol_speed = std::min(vol_speed, max_vol);

    // Clamp injection temp to safety max
    if (injection_temp > 0) {
        int max_temp = config.nozzle_temperature_range_high.get_at(extruder_id);
        if (max_temp > 0)
            injection_temp = std::min(injection_temp, max_temp);
    }

    // Get current printing temperature for restore
    int print_temp = config.nozzle_temperature.get_at(extruder_id);

    // Convert volumetric speed to filament feedrate
    double feedrate_mms = vol_speed / filament_area;
    double feedrate_mmmin = feedrate_mms * 60.0;

    // --- Heat up ---
    // Skip if filament switch already handled temperature via set_extruder()
    bool temp_changed = (injection_temp > 0 && injection_temp != print_temp
                         && !need_filament_switch);
    if (temp_changed) {
        gcode += "; Magma injection: heating\n";
        gcode += gcodegen.retract(false, false);

        if (park_enabled) {
            double park_z = layer_z + 2.0;
            gcode += gcodegen.writer().travel_to_z(park_z, "park z-hop for temp change");
        }

        gcode += gcodegen.writer().set_temperature(injection_temp, true);
    }

    // --- Injection loop ---
    float display_dim = tube_map.interior_width();

    constexpr double slam_depth = 0.1;  // mm, hardcoded for safety

    for (const auto& pt : points) {
        double volume = pt.volume_mm3 * fill_factor;
        if (volume <= 0)
            continue;

        // Travel to injection point (before role tag so unretract classifies correctly)
        Point scaled_pos(scale_(pt.position.x()), scale_(pt.position.y()));
        gcode += gcodegen.travel_to(scaled_pos, erMagmaInjection, "move to injection point");

        // Unretract before injection — happens before the role tag so the
        // GCodeProcessor sees this as a normal Unretract, not an Extrude.
        gcode += gcodegen.unretract();

        // Z-slam: lower nozzle into surface to seal against hole
        if (z_slam) {
            sprintf(buf, "G1 Z%.3f F600 ; z-slam seal\n", layer_z - slam_depth);
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

        // Build tube visualization waypoints (simplified with 3D RDP)
        const auto& pair = tube_map.u_tube_pairs()[pt.pair_index];
        auto waypoints = build_tube_viz_waypoints(
            tube_map, pair, layer_z, pt.window_center_layer);

        // Emit tube visualization metadata for GCodeProcessor
        if (!waypoints.empty())
            gcode += format_tube_viz_comment(waypoints, tube_map.interior_width());

        // Split injection into per-segment G1 E commands so the preview
        // slider shows progressive tube filling.  Each G1 gets its own
        // G-code line number, giving it a separate slider position.
        double filament_length = volume * e_per_mm3;
        gcode += gcodegen.writer().set_speed(feedrate_mmmin);
        Vec2d xy(gcodegen.writer().get_position().x(),
                 gcodegen.writer().get_position().y());

        if (waypoints.size() >= 2) {
            double total_path_len = 0;
            for (size_t i = 1; i < waypoints.size(); ++i)
                total_path_len += (waypoints[i] - waypoints[i - 1]).norm();

            for (size_t i = 1; i < waypoints.size(); ++i) {
                double seg_len = (waypoints[i] - waypoints[i - 1]).norm();
                double seg_e = (total_path_len > 0)
                    ? filament_length * (seg_len / total_path_len)
                    : filament_length / double(waypoints.size() - 1);
                gcode += gcodegen.writer().extrude_to_xy(
                    xy, seg_e, "injection segment");
            }
        } else {
            sprintf(buf, "Magma injection %.2f mm3", volume);
            gcode += gcodegen.writer().extrude_to_xy(
                xy, filament_length, std::string(buf));
        }

        // Dwell for air displacement
        if (dwell_ms > 0) {
            sprintf(buf, "G4 P%d ; injection dwell\n", dwell_ms);
            gcode += buf;
        }

        // Z-slam release: return to normal layer height
        if (z_slam) {
            sprintf(buf, "G1 Z%.3f F600 ; z-slam release\n", layer_z);
            gcode += buf;
        }

        // Retract after injection — relieves ooze pressure, triggers wipe
        gcode += gcodegen.retract(false, false);
    }

    // Reset forced dimensions
    sprintf(buf, ";%s0\n",
            GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Height).c_str());
    gcode += buf;
    sprintf(buf, ";%s0\n",
            GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Width).c_str());
    gcode += buf;

    // --- Ironing of tube ends ---
    if (iron_tube_ends) {
        // Emit ironing role tag
        sprintf(buf, ";%s%s\n",
                GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Role).c_str(),
                ExtrusionEntity::role_to_string(erIroning).c_str());
        gcode += buf;

        // Ironing parameters: use user's settings if ironing is configured,
        // otherwise use sensible defaults
        double ir_flow_pct = 0.10;
        double ir_spacing  = 0.1;
        double ir_speed    = 15.0 * 60.0;  // 15 mm/s → mm/min

        if (config.ironing_type.value != IroningType::NoIroning) {
            ir_flow_pct = config.ironing_flow.value / 100.0;
            ir_spacing  = config.ironing_spacing.value;
            ir_speed    = config.ironing_speed.value * 60.0;  // mm/s → mm/min
        }

        float nozzle_d = config.nozzle_diameter.get_at(extruder_id);
        double ir_height = actual_layer_height * ir_flow_pct;
        double ir_flow_mm3_per_mm = nozzle_d * ir_height;
        double ir_e_per_mm = ir_flow_mm3_per_mm * e_per_mm3;

        // Iron each injection hole with serpentine parallel lines
        for (const auto& pt : points) {
            double radius = tube_map.interior_width() / 2.0 - nozzle_d / 2.0;
            if (radius <= 0.05)
                continue;

            bool left_to_right = true;
            for (double dy = -radius; dy <= radius + 0.001; dy += ir_spacing) {
                double r2 = radius * radius - dy * dy;
                if (r2 <= 0)
                    continue;
                double half_chord = std::sqrt(r2);

                double x0 = pt.position.x() + (left_to_right ? -half_chord : half_chord);
                double x1 = pt.position.x() + (left_to_right ? half_chord : -half_chord);
                double y  = pt.position.y() + dy;

                // Travel to line start
                Point start_s(scale_(x0), scale_(y));
                gcode += gcodegen.travel_to(start_s, erIroning, "iron start");
                gcode += gcodegen.unretract();

                // Extrude ironing line
                double line_len = 2.0 * half_chord;
                double e_val = line_len * ir_e_per_mm;
                gcode += gcodegen.writer().set_speed(ir_speed);
                gcode += gcodegen.writer().extrude_to_xy(
                    Vec2d(x1, y), e_val, "tube iron");

                left_to_right = !left_to_right;
            }
            gcode += gcodegen.retract(false, false);
        }
    }

    // --- Cool down ---
    if (temp_changed) {
        gcode += "; Magma injection: cooling\n";
        gcode += gcodegen.retract(false, false);

        if (park_enabled) {
            double park_z = layer_z + 2.0;
            gcode += gcodegen.writer().travel_to_z(park_z, "park z-hop for temp restore");
        }

        gcode += gcodegen.writer().set_temperature(print_temp, true);
    }

    // --- Switch back to original filament ---
    if (need_filament_switch) {
        gcode += "; Magma injection: restoring print filament\n";
        gcode += gcodegen.set_extruder(original_filament, layer_z);
    }

    return gcode;
}

} // namespace magma
} // namespace Slic3r
