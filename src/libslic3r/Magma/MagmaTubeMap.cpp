#include "MagmaTubeMap.hpp"
#include "MagmaResolved.hpp"
#include "MagmaTubeSolver.hpp"
#include "MagmaPatterns.hpp"

#include "../Layer.hpp"
#include "../Flow.hpp"
#include "../Print.hpp"
#include "../Slicing.hpp"
#include "../ClipperUtils.hpp"
#include "../Polygon.hpp"
#include "../ExPolygon.hpp"
#include "../Surface.hpp"
#include "../ExtrusionEntityCollection.hpp"   // polygons_covered_by_width (measured volume)
#include "../I18N.hpp"
#include "../format.hpp"

#define L(s) Slic3r::I18N::translate(s)

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <set>
#include <boost/log/trivial.hpp>

namespace Slic3r {
namespace magma {

// ============================================================================
// CellPresence
// ============================================================================

bool CellPresence::present(int layer_id) const
{
    if (layer_id < first_layer || layer_id > last_layer)
        return false;
    int idx = layer_id - first_layer;
    return idx >= 0 && idx < int(layers.size()) && layers[idx];
}

double CellPresence::area(int layer_id) const
{
    if (! present(layer_id))
        return 0.0;
    int idx = layer_id - first_layer;
    if (idx < 0 || idx >= int(areas.size()))
        return 0.0;
    return areas[idx];
}

double CellPresence::distance(int layer_id) const
{
    if (! present(layer_id))
        return 0.0;
    int idx = layer_id - first_layer;
    if (idx < 0 || idx >= int(distances.size()))
        return 0.0;
    return distances[idx];
}

double CellPresence::opening_radius(int layer_id) const
{
    if (! present(layer_id))
        return 0.0;
    int idx = layer_id - first_layer;
    if (idx < 0 || idx >= int(opening_radii.size()))
        return 0.0;
    return opening_radii[idx];
}

std::optional<Vec2d> CellPresence::injection_point(int layer_id) const
{
    if (! present(layer_id))
        return std::nullopt;
    int idx = layer_id - first_layer;
    if (idx < 0 || idx >= int(injection_pts.size()))
        return std::nullopt;
    return injection_pts[idx];
}

void CellPresence::mark_present(int layer_id, double area_val, double dist_val,
                               double opening_r, const Vec2d &inj_pt)
{
    if (first_layer == INT_MAX) {
        // First time: initialize
        first_layer = layer_id;
        last_layer  = layer_id;
        layers.assign(1, true);
        areas.assign(1, area_val);
        distances.assign(1, dist_val);
        opening_radii.assign(1, opening_r);
        injection_pts.assign(1, inj_pt);
        return;
    }

    // Expand range if needed
    if (layer_id < first_layer) {
        int extend = first_layer - layer_id;
        layers.insert(layers.begin(), extend, false);
        areas.insert(areas.begin(), extend, 0.0);
        distances.insert(distances.begin(), extend, 0.0);
        opening_radii.insert(opening_radii.begin(), extend, 0.0);
        injection_pts.insert(injection_pts.begin(), extend, Vec2d(0.0, 0.0));
        first_layer = layer_id;
    }
    if (layer_id > last_layer) {
        int extend = layer_id - last_layer;
        layers.insert(layers.end(), extend, false);
        areas.insert(areas.end(), extend, 0.0);
        distances.insert(distances.end(), extend, 0.0);
        opening_radii.insert(opening_radii.end(), extend, 0.0);
        injection_pts.insert(injection_pts.end(), extend, Vec2d(0.0, 0.0));
        last_layer = layer_id;
    }

    int idx = layer_id - first_layer;
    layers[idx]        = true;
    areas[idx]         = area_val;
    distances[idx]     = dist_val;
    opening_radii[idx] = opening_r;
    injection_pts[idx] = inj_pt;
}

// ============================================================================
// MagmaTubeMap — Statistics
// ============================================================================

int MagmaTubeMap::num_solid_cells() const
{
    int count = 0;
    for (const auto &kv : m_cell_pair_index)
        if (kv.second.empty()) ++count;
    return count;
}

// ============================================================================
// MagmaTubeMap — Per-layer data accessors and helpers
// ============================================================================

double MagmaTubeMap::print_z(int layer_id) const
{
    if (valid_layer(layer_id))
        return m_layer_data[layer_id].print_z;
    return 0.0;
}

double MagmaTubeMap::layer_height_at(int layer_id) const
{
    if (valid_layer(layer_id))
        return m_layer_data[layer_id].height;
    return double(m_layer_height);
}

int MagmaTubeMap::window_center_layer(const UTubePair& pair) const
{
    double start_bottom = print_z(pair.pair_start_layer)
                        - layer_height_at(pair.pair_start_layer);
    double window_mid_z = (start_bottom + pair.window_end_z) / 2.0;
    for (int l = pair.pair_start_layer; l <= pair.pair_end_layer; ++l) {
        if (print_z(l) >= window_mid_z)
            return l;
    }
    return pair.pair_start_layer;
}

double MagmaTubeMap::span_height_mm(int start_layer, int end_layer) const
{
    if (end_layer < start_layer || ! valid_layer(start_layer) || ! valid_layer(end_layer))
        return 0.0;
    // print_z is cumulative (top of layer), so span height =
    // top of end_layer minus bottom of start_layer.
    return m_layer_data[end_layer].print_z - m_layer_data[start_layer].bottom_z();
}

// ============================================================================
// MagmaTubeMap — Build
// ============================================================================

std::unique_ptr<MagmaTubeMap> MagmaTubeMap::build(
    const std::vector<Layer*> &layers,
    const PrintRegionConfig &config,
    const PrintConfig &print_config,
    const PrintObjectConfig &obj_config,
    const SlicingParameters &slicing_params,
    ProgressFn progress_fn,
    ThrowIfCanceled throw_if_canceled)
{
    auto t_start = std::chrono::high_resolution_clock::now();

    auto map = std::unique_ptr<MagmaTubeMap>(new MagmaTubeMap());

    MagmaResolved res;
    if (! resolve_magma(config, obj_config, print_config, res)) {
        // Unreachable from PrintObject, which only calls this for a Magma region; log it here
        // rather than leave only FillMagmaBase's later "no tube map" error.
        BOOST_LOG_TRIVIAL(error)
            << "MagmaTubeMap::build called for a region whose effective pattern is not a Magma "
               "pattern; no tube map will be built.";
        return nullptr;
    }

    map->m_injection_extruder = res.injection_extruder;
    map->m_line_width     = static_cast<float>(res.line_width);
    map->m_pattern        = res.pattern;
    map->m_geometry       = res.geometry;
    map->m_interior_width = static_cast<float>(res.interior_width);
    map->m_cell_spacing   = res.cell_spacing;

    // Nominal layer height, the fallback for layer ids outside the per-layer table.
    map->m_layer_height = static_cast<float>(obj_config.layer_height.value);
    map->m_solver_mode    = obj_config.magma_tube_solver_mode.value;
    map->m_solver_timeout = obj_config.magma_solver_timeout.value;

    // Per-layer height/Z tables (adaptive layer heights). Must precede the spiral params,
    // which use m_min_layer_height.
    {
        // m_layer_data is indexed by absolute Layer::id(), which starts after the raft layers.
        // Rows below m_first_layer_id stay zero-filled and must not be read as real layers.
        int max_layer_id = 0;
        int min_layer_id = INT_MAX;
        for (const Layer *l : layers) {
            max_layer_id = std::max(max_layer_id, static_cast<int>(l->id()));
            min_layer_id = std::min(min_layer_id, static_cast<int>(l->id()));
        }
        map->m_first_layer_id = (min_layer_id == INT_MAX) ? 0 : min_layer_id;
        map->m_num_layers = max_layer_id + 1;
        map->m_layer_data.resize(map->m_num_layers);
        for (const Layer *l : layers) {
            int lid = static_cast<int>(l->id());
            map->m_layer_data[lid].print_z = l->print_z;
            map->m_layer_data[lid].height  = l->height;
        }
        // Smallest layer height this object actually prints -- not
        // slicing_params.min_layer_height, which is the machine floor.
        {
            const double min_h = min_positive_layer_height(map->m_layer_data,
                                                           map->m_first_layer_id);
            map->m_min_layer_height = static_cast<float>(min_h > 0.0 ? min_h : map->m_layer_height);
        }
    }

    // The helix angle is capped at the thinnest layer.
    bool spiral_enabled = obj_config.magma_spiral_interlock.value;
    map->m_spiral_params = compute_spiral_params(map->m_interior_width, map->m_line_width,
                                                  map->m_min_layer_height, spiral_enabled);

    // Window spec (handles window height auto-calculation).
    map->m_window_spec = WindowSpec::from_config(
        *map->m_geometry,
        static_cast<float>(obj_config.magma_window_height_mm.value),
        map->m_interior_width,
        map->m_line_width,
        map->m_layer_height
    );

    // Min tube height: the tube must clear its own window by a full window height, so the
    // shared wall above it keeps the pair two legs rather than one chamber.
    map->m_min_tube_height_mm = magma_min_tube_height(map->m_window_spec.window_height_mm);

    // Max tube height from config, floored at exactly m_min_tube_height_mm with no headroom:
    // it governs injection duration, and whether a layer boundary lands inside [min, max] is
    // the solver's job.
    {
        double user_tube_mm = std::max(1.0, obj_config.magma_tube_height.value);
        map->m_max_tube_height_mm = std::max(user_tube_mm, map->m_min_tube_height_mm);

        if (user_tube_mm < map->m_min_tube_height_mm) {
            map->m_warning_messages.push_back(Slic3r::format(
                "Magma max tube height (%.1fmm) is too short to fit a complete U-tube pair. "
                "Minimum is %.1fmm (twice the window height). "
                "Using that minimum.",
                user_tube_mm, map->m_min_tube_height_mm));
        }
    }

    // Injection speed vs filament max volumetric speed.
    {
        double inj_speed = obj_config.magma_injection_speed.value;
        double max_vol = print_config.filament_max_volumetric_speed.get_at(
            (unsigned int)res.injection_extruder);
        if (max_vol > 0 && inj_speed > max_vol) {
            map->m_warning_messages.push_back(Slic3r::format(
                "Magma injection speed (%.1f mm\u00b3/s) exceeds the filament's "
                "max volumetric speed (%.1f mm\u00b3/s). It will be capped automatically.",
                inj_speed, max_vol));
        }
    }

    // Per-layer lattice cache, each with that layer's spiral offset.
    for (int i = map->m_first_layer_id; i < map->m_num_layers; ++i)
        map->m_layer_data[i].lattice = lattice_for_layer(map->m_pattern, map->m_cell_spacing, map->m_spiral_params, i, map->m_line_width);

    map->m_injection_edge_pref = obj_config.magma_injection_edge_pref.value;

    // Cap-seal clearance gate for scan_layers: the injection nozzle flat seats at the injection
    // point, so a tube layer needs >= flat/2 clearance from it to the nearest boundary. Exactly
    // flat/2, with no margin, so boundary cells just over the limit are kept.
    map->m_min_cap_clearance = 0.5 * res.injection_nozzle_flat;

    map->scan_layers(layers);
    map->assign_tubes(progress_fn, throw_if_canceled);
    map->assign_extra_vents();   // tri-hex: add manifold legs (no-op for other patterns)
    map->precompute_window_end_z();
    map->precompute_injection_data();

    // pair.volume_mm3 and the cap-layer list stay empty until measure_volumes() runs after
    // PrintObject::infill().

    // Magma is enabled but produced no tubes: warn with the most likely reason.
    if (map->m_pairs.empty()) {
        // Tallest contiguous single-cell column. If even this is below the minimum tube
        // height, no tube can form regardless of pairing; otherwise the failure is in
        // pairing (adjacency / shared height), not part height.
        double tallest_run_mm = 0.0;
        for (const auto &kv : map->m_cells) {
            const CellPresence &cp = kv.second;
            if (cp.first_layer == INT_MAX) continue;
            int run_start = -1;
            for (int L = cp.first_layer; L <= cp.last_layer; ++L) {
                bool present = cp.present(L);
                if (present && run_start < 0) run_start = L;
                if ((!present || L == cp.last_layer) && run_start >= 0) {
                    int run_end = present ? L : L - 1;
                    tallest_run_mm = std::max(tallest_run_mm,
                                              map->span_height_mm(run_start, run_end));
                    run_start = -1;
                }
            }
        }

        std::string reason;
        if (map->m_cells.empty()) {
            reason =
                "no Magma cells fit inside the part — the cell size is too large for this "
                "object's cross-section. Reduce the Magma interior width or infill spacing.";
        } else if (tallest_run_mm + 1e-6 < map->m_min_tube_height_mm) {
            reason = Slic3r::format(
                "the reinforced region is only %.1fmm tall, below the %.1fmm minimum for one "
                "U-tube (twice the window height). Use a taller part, or reduce the "
                "Magma interior width to shrink the window (and the minimum).",
                tallest_run_mm, map->m_min_tube_height_mm);
        } else {
            reason =
                "cells are tall enough but none could be paired into U-tubes — neighbouring "
                "cells may not share enough height or area. Try a smaller cell size, denser "
                "infill spacing, or a more uniform part shape.";
        }
        map->m_warning_messages.push_back("Magma produced no reinforcement tubes: " + reason);
    }

    // Tri-hex: warn when the triangle vents are too small to stay open. A vent is an
    // equilateral triangle of side cell_spacing/sqrt3, so its open inscribed diameter is
    // cell_spacing/3 - line_width (must match cell_bore_at()). Below ~1mm, line overlap and
    // ordinary print imperfections seal it.
    if (map->m_pattern == ipMagmaTriHex) {
        const double vent_dia = map->m_cell_spacing / 3.0 - double(map->m_line_width);
        if (vent_dia < 1.0) {
            map->m_warning_messages.push_back(Slic3r::format(
                "Magma Tri-hex vents are only %.2f mm across after line overlap — they may seal "
                "off and block injection. Tri-hex needs large tubes: use an interior width of "
                "about 4 mm or more (or switch to Magma Honeycomb / Rectilinear) for reliable vents.",
                vent_dia));
        }
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

    BOOST_LOG_TRIVIAL(info) << "MagmaTubeMap built in " << ms << "ms: "
        << map->num_cells() << " cells, "
        << map->num_pairs() << " pairs, "
        << map->num_solid_cells() << " solid | "
        << "max_tube_height=" << map->m_max_tube_height_mm << "mm"
        << " min_tube_height=" << map->m_min_tube_height_mm << "mm";

    return map;
}

// ============================================================================
// MagmaTubeMap — tube_opening_diameter
// ============================================================================

double MagmaTubeMap::tube_opening_diameter() const
{
    return m_geometry->opening_diameter(m_cell_spacing, m_line_width);
}

double MagmaTubeMap::cap_opening_diameter(const UTubePair &pair) const
{
    // Prefer the deposited cap cavity measured post-fill: the scan's cell_inset model
    // over-reads hex openings, whose walls are doubled.
    if (pair.measured_cap_opening_dia > 0.0)
        return pair.measured_cap_opening_dia;
    auto it = m_cells.find(pair.cell_a);
    if (it != m_cells.end()) {
        double r = it->second.opening_radius(pair.pair_end_layer);
        if (r > 0.0)
            return 2.0 * r;
    }
    return tube_opening_diameter();
}

double MagmaTubeMap::cell_bore_at(const CellId &cell, int layer) const
{
    // Unclipped bore: the interior width for kind 0 (every single-kind pattern, and the
    // tri-hex hub); for a tri-hex vent, the inscribed diameter of the inset triangle,
    // spacing/3 - lw, which must match the vent-width warning in build().
    double ideal_bore;
    if (cell.kind == 0) {
        ideal_bore = m_interior_width;
    } else {
        double e     = m_cell_spacing * INV_SQRT3;
        double inset = e - double(m_line_width) * SQRT3;
        ideal_bore = inset > 0.0 ? inset * INV_SQRT3 : m_interior_width;
    }
    // Area scales with bore², hence the sqrt. max_area approximates the unclipped area.
    auto it = m_cells.find(cell);
    if (it == m_cells.end())
        return ideal_bore;
    const CellPresence &p = it->second;
    double a = p.area(layer);
    if (a <= 0.0)
        return ideal_bore;
    double amax = 0.0;
    for (double av : p.areas)
        amax = std::max(amax, av);
    if (amax <= 0.0)
        return ideal_bore;
    return ideal_bore * std::sqrt(std::min(1.0, a / amax));
}

const MagmaLattice& MagmaTubeMap::topology_lattice() const
{
    if (!m_topology_lattice)
        m_topology_lattice = make_magma_lattice(m_pattern, m_cell_spacing, 0.0, 0.0, m_line_width);
    return *m_topology_lattice;
}

// ============================================================================
// cell_inset_polygons — a cell's hollow open cross-section
// ============================================================================
//
// The cell's wall outline (cell_corners) inset by half the bead width. offset_ex can split
// or empty a thin cell, hence ExPolygons.
static ExPolygons cell_inset_polygons(const MagmaLattice &lattice,
                                      const CellId &cell,
                                      double half_line_width)
{
    std::vector<Vec2d> corners = lattice.cell_corners(cell);
    Polygon poly;
    poly.points.reserve(corners.size());
    for (const Vec2d &c : corners)
        poly.points.emplace_back(scale_(c.x()), scale_(c.y()));
    return offset_ex(poly, -scale_(half_line_width));
}

// ============================================================================
// MagmaTubeMap — scan_layers
// ============================================================================

void MagmaTubeMap::scan_layers(const std::vector<Layer*> &layers)
{
    const double half_line_width = m_line_width * 0.5;

    // A cell is present on a layer if its clipped open area is >= 70% of its kind's ideal.
    // Deliberately loose: injection volume is measured from the actual clipped area, and
    // reliability depends on the cap-layer seal, not per-layer fullness. It also keeps boundary
    // cells on the first layer, whose fill region is shrunk by the wider first-layer wall.
    const double presence_area_frac = 0.70;

    // A cell whose centre is this far inside the zone has its whole inscribed bore
    // (diameter = interior_width) inside it.
    const coord_t interior_inset = scale_(m_interior_width * 0.5);

    // Cell identity comes from the zero-offset lattice so it is stable across layers; the
    // spiral-offset per-layer lattice is used for position checks.
    const MagmaLattice &ref_lattice = topology_lattice();

    // Per-kind ideal {open area (scaled^2), opening radius (mm)}, cached on first use and
    // measured with the same cell_inset_polygons as boundary cells. Per kind so tri-hex vents
    // are not sized as hubs. Opening radius = farthest inset corner from the cell centre.
    std::pair<double, double> kind_ideal[3];
    bool kind_done[3] = { false, false, false };
    auto ideal_for = [&](const CellId &cell) -> std::pair<double, double> {
        const uint8_t k = cell.kind;
        if (k < 3 && kind_done[k]) return kind_ideal[k];
        ExPolygons inset = cell_inset_polygons(ref_lattice, cell, half_line_width);
        double area = 0.0;
        for (const ExPolygon &ep : inset)
            area += std::abs(ep.area());
        Vec2d  ctr = ref_lattice.cell_center(cell);
        Point  ctr_pt(scale_(ctr.x()), scale_(ctr.y()));
        double opr_scaled = 0.0;
        for (const ExPolygon &ep : inset)
            for (const Point &p : ep.contour.points)
                opr_scaled = std::max(opr_scaled, (ctr_pt - p).cast<double>().norm());
        std::pair<double, double> val{ area, unscale<double>(opr_scaled) };
        if (k < 3) { kind_ideal[k] = val; kind_done[k] = true; }
        return val;
    };

    for (int i = 0; i < int(layers.size()); ++i) {
        const Layer *layer = layers[i];
        // Layer::id(), not the array index: it includes the raft offset, and FillMagma and
        // GCode query with it.
        const int layer_id = static_cast<int>(layer->id());

        const MagmaLattice &layer_lattice = *m_layer_data[layer_id].lattice;

        // Collect the surfaces where the lattice is printed.
        ExPolygons zone_regions;
        for (const LayerRegion *layerm : layer->regions()) {
            const PrintRegionConfig &region_config = layerm->region().config();
            for (const Surface &surface : layerm->fill_surfaces.surfaces) {
                // Non-dual Magma: all sparse infill is lattice. Dual infill: only layers that
                // have a zone boundary (has_dual_infill_zone()).
                if (surface.surface_type == stInternal
                           && (layerm->has_dual_infill_zone() || is_magma_pattern(region_config.sparse_infill_pattern.value))) {
                    zone_regions.push_back(surface.expolygon);
                }
            }
        }

        zone_regions = union_ex(zone_regions);
        if (zone_regions.empty())
            continue;

        // Cells centred in here need no per-cell clipping.
        ExPolygons interior_region = offset_ex(zone_regions, -interior_inset);

        BoundingBox bbox = get_extents(zone_regions);
        // Margin for the spiral offset moving cells.
        bbox.offset(scale_(m_interior_width));
        std::vector<CellId> cells = ref_lattice.enumerate_cells(bbox);

        for (const CellId &cell : cells) {
            Vec2d center_mm = layer_lattice.cell_center(cell);
            Point center_pt(scale_(center_mm.x()), scale_(center_mm.y()));

            // Fast path: the whole bore is inside the zone.
            bool is_interior = false;
            for (const ExPolygon &ep : interior_region) {
                if (ep.contains(center_pt)) {
                    is_interior = true;
                    break;
                }
            }

            if (is_interior) {
                Point proj = projection_onto(zone_regions, center_pt);
                double dist_mm = unscale<double>((center_pt - proj).cast<double>().norm());
                // Same clearance gate as the boundary path; catches an interior width too
                // small for even an unclipped cell to seat the nozzle.
                if (dist_mm < m_min_cap_clearance)
                    continue;
                // Unclipped: cavity == full cell, so the injection point is the cell centre.
                std::pair<double, double> ideal = ideal_for(cell);
                m_cells[cell].mark_present(layer_id, ideal.first, dist_mm,
                                           ideal.second, center_mm);
                continue;
            }

            bool in_region = false;
            for (const ExPolygon &ep : zone_regions) {
                if (ep.contains(center_pt)) {
                    in_region = true;
                    break;
                }
            }

            if (!in_region)
                continue;

            // Boundary cell: open area = the inset cell clipped to the zone.
            ExPolygons inset = cell_inset_polygons(layer_lattice, cell, half_line_width);
            if (inset.empty())
                continue;

            ExPolygons clipped = intersection_ex(inset, zone_regions);
            if (clipped.empty())
                continue;

            double area = 0.0;
            for (const ExPolygon &ep : clipped)
                area += std::abs(ep.area());

            if (area < presence_area_frac * ideal_for(cell).first)
                continue;

            // Injection point: centroid of the largest clipped piece, which moves away from the
            // part wall for more nozzle clearance. Falls back to the cell centre if a concave
            // clip puts the centroid outside the opening.
            const ExPolygon *piece = &clipped.front();
            for (const ExPolygon &ep : clipped)
                if (std::abs(ep.area()) > std::abs(piece->area())) piece = &ep;
            Point inj_pt = piece->contour.centroid();
            if (!piece->contains(inj_pt))
                inj_pt = center_pt;
            Vec2d inj_mm(unscale<double>(inj_pt.x()), unscale<double>(inj_pt.y()));

            // Clearance gate (see m_min_cap_clearance), measured at the injection point itself:
            // a spike toward it fails the gate even when the area passes.
            double clearance_mm = unscale<double>(
                (inj_pt - projection_onto(clipped, inj_pt)).cast<double>().norm());
            if (clearance_mm < m_min_cap_clearance)
                continue;

            // Opening radius: farthest point of the clipped opening from the injection point,
            // which the seal must cover (see cap_opening_diameter()).
            double opening_r_scaled = 0.0;
            for (const ExPolygon &ep : clipped)
                for (const Point &p : ep.contour.points)
                    opening_r_scaled = std::max(opening_r_scaled,
                                                (inj_pt - p).cast<double>().norm());
            double opening_r = unscale<double>(opening_r_scaled);

            // Cell centre -> nearest boundary, for the injection-side preference.
            Point proj = projection_onto(zone_regions, center_pt);
            double dist_mm = unscale<double>((center_pt - proj).cast<double>().norm());

            m_cells[cell].mark_present(layer_id, area, dist_mm, opening_r, inj_mm);
        }
    }

    // Log how far cells reach toward the object's top layer.
    {
        int max_layer_id = layers.empty() ? -1
            : static_cast<int>(layers.back()->id());
        int cells_reaching_top = 0;
        int min_last = INT_MAX, max_last = -1;
        for (const auto &[cell, presence] : m_cells) {
            if (presence.last_layer == max_layer_id) ++cells_reaching_top;
            min_last = std::min(min_last, presence.last_layer);
            max_last = std::max(max_last, presence.last_layer);
        }
        BOOST_LOG_TRIVIAL(info) << "MagmaScanSpan top_layer_id=" << max_layer_id
            << " cells=" << m_cells.size()
            << " reaching_top=" << cells_reaching_top
            << " last_layer[min=" << (min_last == INT_MAX ? -1 : min_last)
            << " max=" << max_last << "]";
    }
}

// ============================================================================
// MagmaTubeMap — assign_tubes (delegates to CP-SAT solver)
// ============================================================================

void MagmaTubeMap::assign_tubes(ProgressFn progress_fn, ThrowIfCanceled throw_if_canceled)
{
    MagmaTubeSolver solver(topology_lattice(), m_cells, m_layer_data, m_first_layer_id,
                           m_min_tube_height_mm, m_max_tube_height_mm,
                           m_num_layers,
                           m_solver_mode, m_solver_timeout);
    solver.solve(m_pairs, m_cell_pair_index, progress_fn, throw_if_canceled);

    // A failed block may be infeasible or timed out; the solver does not report which, so the
    // message does not claim a cause.
    if (solver.failed_block_count() > 0) {
        m_warning_messages.push_back(Slic3r::format(
            "Magma tube solver: %1% block(s) could not be solved. The previous result was "
            "kept for those blocks, so reinforcement there may be less than intended. "
            "Raising the solver timeout may help if the model is merely slow.",
            solver.failed_block_count()));
    }
}

// ============================================================================
// MagmaTubeMap — assign_extra_vents (tri-hex manifold legs)
// ============================================================================
//
// The solver gives each hub-tube one primary vent. This additive pass gives triangle vents
// extra legs on bordering hub-tubes: a hub-tube is a candidate only if the vent is present
// and not a primary on every layer of its range (a gap would trap air). Each vent takes a
// max-coverage non-overlapping subset of its candidates by weighted interval scheduling.
// Hubs accept any number of legs, and the solver already made every hub feasible, so this
// cannot strand a hub.

void MagmaTubeMap::assign_extra_vents()
{
    if (m_pattern != ipMagmaTriHex || m_pairs.empty())
        return;

    const MagmaLattice &lat = topology_lattice();

    auto hub_of  = [](const UTubePair &p) { return p.cell_a.kind == THK_HEX ? p.cell_a : p.cell_b; };
    auto vent_of = [](const UTubePair &p) { return p.cell_a.kind == THK_HEX ? p.cell_b : p.cell_a; };

    // Index hub-tubes by their hub cell; collect each vent's primary-claimed ranges.
    std::unordered_map<CellId, std::vector<int>, CellIdHash> tubes_by_hub;
    std::unordered_map<CellId, std::vector<std::pair<int, int>>, CellIdHash> primary_claim;
    for (int pi = 0; pi < int(m_pairs.size()); ++pi) {
        const UTubePair &p = m_pairs[pi];
        tubes_by_hub[hub_of(p)].push_back(pi);
        primary_claim[vent_of(p)].push_back({ p.pair_start_layer, p.pair_end_layer });
    }

    int extra_legs_added = 0;

    for (const auto &cell_kv : m_cells) {
        const CellId &V = cell_kv.first;
        if (V.kind == THK_HEX)              // hubs are not vents
            continue;
        const CellPresence &pres = cell_kv.second;

        auto it_claim = primary_claim.find(V);
        auto available_at = [&](int L) -> bool {
            if (!pres.present(L))
                return false;
            if (it_claim != primary_claim.end())
                for (const auto &rng : it_claim->second)
                    if (L >= rng.first && L <= rng.second)
                        return false;       // busy as this layer's primary
            return true;
        };

        // Feasible candidate hub-tubes on V's bordering hubs (excluding tubes where V is
        // already the primary) whose ENTIRE range is available.
        struct Cand { int pi; int start; int cap; };
        std::vector<Cand> cands;
        for (const CellId &H : lat.neighbors(V)) {
            auto it = tubes_by_hub.find(H);
            if (it == tubes_by_hub.end())
                continue;
            for (int pi : it->second) {
                const UTubePair &p = m_pairs[pi];
                if (vent_of(p) == V)
                    continue;
                bool ok = true;
                for (int L = p.pair_start_layer; L <= p.pair_end_layer; ++L)
                    if (!available_at(L)) { ok = false; break; }
                if (ok)
                    cands.push_back({ pi, p.pair_start_layer, p.pair_end_layer });
            }
        }
        if (cands.empty())
            continue;

        // Weighted interval scheduling: max total layers covered by a non-overlapping
        // subset. Sort by cap, DP with latest-compatible-predecessor, then backtrack.
        std::sort(cands.begin(), cands.end(),
                  [](const Cand &a, const Cand &b) { return a.cap < b.cap; });
        const int n = int(cands.size());
        auto wt = [&](int i) { return cands[i].cap - cands[i].start + 1; };  // layers covered
        std::vector<int> pred(n, -1);
        for (int i = 0; i < n; ++i)
            for (int j = i - 1; j >= 0; --j)
                if (cands[j].cap < cands[i].start) { pred[i] = j; break; }
        std::vector<int> best(n + 1, 0);   // best[i] over cands[0..i-1]
        for (int i = 1; i <= n; ++i) {
            int incl = wt(i - 1) + (pred[i - 1] >= 0 ? best[pred[i - 1] + 1] : 0);
            best[i] = std::max(best[i - 1], incl);
        }
        for (int i = n; i > 0; ) {
            int incl = wt(i - 1) + (pred[i - 1] >= 0 ? best[pred[i - 1] + 1] : 0);
            if (incl > best[i - 1]) {       // candidate i-1 is in the optimal set
                const int pi = cands[i - 1].pi;
                m_pairs[pi].extra_vents.push_back(V);
                // Register the leg like the solver registers cell_a/cell_b.
                m_cell_pair_index[V].push_back(pi);
                ++extra_legs_added;
                i = (pred[i - 1] >= 0) ? pred[i - 1] + 1 : 0;
            } else {
                --i;
            }
        }
    }

    BOOST_LOG_TRIVIAL(info) << "MagmaTriHex extra-vent sweep: added " << extra_legs_added
        << " extra vent legs across " << m_pairs.size() << " hub-tubes";
}

// ============================================================================
// MagmaTubeMap — precompute_window_end_z
// ============================================================================

void MagmaTubeMap::precompute_window_end_z()
{
    // In mm, not layers, so it holds under variable layer height. A layer is in the window
    // if its bottom is below window_end_z.
    const double wh_mm = m_window_spec.window_height_mm;

    for (UTubePair &pair : m_pairs) {
        double start_bottom = m_layer_data[pair.pair_start_layer].bottom_z();
        pair.window_end_z = start_bottom + wh_mm;
    }
}

// ============================================================================
// MagmaTubeMap — precompute_injection_data
// ============================================================================

void MagmaTubeMap::precompute_injection_data()
{
    for (UTubePair &pair : m_pairs) {
        int cap_layer = pair.pair_end_layer;

        // Choose the injection (cap) side.
        if (pair.cell_a.kind != pair.cell_b.kind) {
            // Tri-hex: always inject into the hexagon hub; edge preference does not apply.
            if (pair.cell_a.kind != THK_HEX)
                std::swap(pair.cell_a, pair.cell_b);
        } else {
            // Single-kind patterns: Interior injects into the cell further from the part edge,
            // Exterior into the nearer one.
            auto it_a = m_cells.find(pair.cell_a);
            auto it_b = m_cells.find(pair.cell_b);
            if (it_a != m_cells.end() && it_b != m_cells.end()) {
                double dist_a = it_a->second.distance(cap_layer);
                double dist_b = it_b->second.distance(cap_layer);
                bool swap = (m_injection_edge_pref == MagmaInjectionEdgePref::Interior)
                    ? (dist_a < dist_b)   // A is closer to edge, want interior → swap
                    : (dist_a > dist_b);  // A is further from edge, want exterior → swap
                if (swap)
                    std::swap(pair.cell_a, pair.cell_b);
            }
        }

        // Injection point from the scan at the cap layer. A paired cell is always present
        // there; if not, warn and use the lattice centre.
        auto it_inj = m_cells.find(pair.cell_a);
        std::optional<Vec2d> inj_pt;
        if (it_inj != m_cells.end())
            inj_pt = it_inj->second.injection_point(cap_layer);
        if (! inj_pt) {
            BOOST_LOG_TRIVIAL(warning)
                << "Magma: no measured injection centroid for cell ("
                << pair.cell_a.a << "," << pair.cell_a.b << "," << pair.cell_a.c
                << ") at its cap layer " << cap_layer
                << "; using the lattice centre.";
        }
        pair.injection_center = inj_pt
            ? *inj_pt
            : m_layer_data[cap_layer].lattice->cell_center(pair.cell_a);
        pair.window_center_layer = window_center_layer(pair);
    }
}

// ============================================================================
// MagmaTubeMap — measure_volumes (from the actual deposited toolpath, post-fill)
// ============================================================================
//
// Per layer, a pair's cavity is (cell_a ∪ cell_b ∪ vents) ∩ zone minus the deposited walls;
// the volume is its area times layer height summed over the tube, less the vertex-overlap
// excess. This accounts for doubled walls, the window gap and part-edge clipping without
// per-pattern corrections. Needs the fills from PrintObject::infill().
void MagmaTubeMap::measure_volumes(const std::vector<Layer*> &layers)
{
    std::unordered_map<int, const Layer*> by_id;
    for (const Layer *L : layers) by_id[int(L->id())] = L;

    // Layer ids touched by any pair.
    std::set<int> used;
    for (const UTubePair &p : m_pairs)
        for (int L = p.pair_start_layer; L <= p.pair_end_layer; ++L) used.insert(L);

    // Per-layer cache (shared by all pairs on a layer): the magma zone + deposited walls.
    std::unordered_map<int, ExPolygons> zone_by_layer, walls_by_layer;
    std::unordered_map<int, double> h_by_layer;
    for (int lid : used) {
        auto it = by_id.find(lid);
        if (it == by_id.end()) continue;
        const Layer *layer = it->second;
        h_by_layer[lid] = double(layer->height);
        ExPolygons zone; Polygons walls;
        for (const LayerRegion *lr : layer->regions()) {
            const PrintRegionConfig &rc = lr->region().config();
            if (!(lr->has_dual_infill_zone() || is_magma_pattern(rc.sparse_infill_pattern.value)))
                continue;
            for (const Surface &s : lr->fill_surfaces.surfaces)
                if (s.surface_type == stInternal) zone.push_back(s.expolygon);
            append(walls, lr->fills.polygons_covered_by_width(0.f));
        }
        zone_by_layer[lid]  = union_ex(zone);
        walls_by_layer[lid] = union_ex(walls);
    }

    // Where infill lines cross, the second bead's material bulges into the void, which the
    // union in polygons_covered_by_width does not capture, so it is subtracted per layer.
    const double excess_unit = m_geometry->vertex_overlap_excess_area(m_line_width);

    m_injection_layer_ids.clear();
    int    zero   = 0;
    double t_meas = 0.0;
    for (UTubePair &pair : m_pairs) {
        // Per-pair overlap excess area (mm^2 per unit height), apportioned per cell. Tri-hex
        // charges by corner count (hub 6, vent 3, /4 into each incident cell); other patterns
        // use the per-cell value (0 for honeycomb, whose lines meet at degree-3 vertices).
        double excess_rate;
        if (m_pattern == ipMagmaTriHex) {
            auto charge = [&](const CellId &c) {
                return 0.25 * ((c.kind == THK_HEX) ? 6.0 : 3.0) * excess_unit;
            };
            excess_rate = charge(pair.cell_a) + charge(pair.cell_b);
            for (const CellId &ev : pair.extra_vents) excess_rate += charge(ev);
        } else {
            excess_rate = excess_unit * double(2 + int(pair.extra_vents.size()));
        }
        double vol = 0.0;
        for (int lid = pair.pair_start_layer; lid <= pair.pair_end_layer; ++lid) {
            auto zit = zone_by_layer.find(lid);
            if (zit == zone_by_layer.end() || zit->second.empty()) continue;
            const MagmaLattice &lat   = lattice_at(lid);
            const ExPolygons   &walls = walls_by_layer[lid];
            // Open cross-section of cells at this layer: corner polygons clipped to the zone,
            // minus the deposited walls.
            auto cell_cavity = [&](const std::vector<CellId> &cs) -> ExPolygons {
                Polygons polys;
                for (const CellId &c : cs) {
                    std::vector<Vec2d> cor = lat.cell_corners(c);
                    if (cor.size() < 3) continue;
                    Polygon poly;
                    for (const Vec2d &v : cor) poly.points.push_back(Point(scale_(v.x()), scale_(v.y())));
                    polys.push_back(std::move(poly));
                }
                ExPolygons z = intersection_ex(union_ex(polys), zit->second);
                return (walls.empty() || z.empty()) ? z : diff_ex(z, walls);
            };

            std::vector<CellId> all_cells{ pair.cell_a, pair.cell_b };
            for (const CellId &ev : pair.extra_vents) all_cells.push_back(ev);
            ExPolygons cavity = cell_cavity(all_cells);
            if (cavity.empty()) continue;
            double a2 = 0.0;
            for (const ExPolygon &ep : cavity) a2 += std::abs(ep.area());
            vol += unscale<double>(unscale<double>(a2)) * h_by_layer[lid];
            vol -= excess_rate * h_by_layer[lid];

            // At the cap layer, the injection cell's cavity gives both the injection center
            // (its centroid) and the seal opening.
            if (lid == pair.pair_end_layer) {
                ExPolygons cap = cell_cavity({ pair.cell_a });
                if (!cap.empty()) {
                    const ExPolygon *big = &cap.front();
                    for (const ExPolygon &ep : cap)
                        if (std::abs(ep.area()) > std::abs(big->area())) big = &ep;
                    Point c = big->contour.centroid();
                    if (!big->contains(c)) {                       // concave clip → cell center
                        Vec2d cc = lat.cell_center(pair.cell_a);
                        c = Point(scale_(cc.x()), scale_(cc.y()));
                    }
                    pair.injection_center = Vec2d(unscale<double>(c.x()), unscale<double>(c.y()));
                    // Seal opening: 2x the farthest cavity point from the injection center,
                    // including wall-gap slivers. Covering that radius seals the whole opening.
                    double r = 0.0;
                    for (const Point &p : big->contour.points)
                        r = std::max(r, (c - p).cast<double>().norm());
                    pair.measured_cap_opening_dia = 2.0 * unscale<double>(r);
                }
            }
        }
        if (vol < 0.0) vol = 0.0;   // overlap can exceed a tiny cell's cavity → not injectable
        t_meas += vol;
        pair.volume_mm3 = vol;
        if (vol > 0.0) m_injection_layer_ids.push_back(pair.pair_end_layer);
        else ++zero;
    }
    sort_remove_duplicates(m_injection_layer_ids);
    BOOST_LOG_TRIVIAL(info) << "Magma measured volumes: " << m_pairs.size() << " pairs, "
        << zero << " zero-volume | total=" << t_meas << "mm3 (overlap@lw="
        << m_line_width << "mm)";
}

// ============================================================================
// MagmaTubeMap — Query Interface
// ============================================================================

bool MagmaTubeMap::window_open_at(const UTubePair &pair, int layer_id) const
{
    // In the pair's layer range and the layer's bottom is below window_end_z.
    return layer_id >= pair.pair_start_layer && layer_id <= pair.pair_end_layer
        && (print_z(layer_id) - layer_height_at(layer_id)) < pair.window_end_z;
}

bool MagmaTubeMap::is_paired(const CellId &cell) const
{
    auto it = m_cell_pair_index.find(cell);
    return (it != m_cell_pair_index.end() && !it->second.empty());
}


} // namespace magma
} // namespace Slic3r
