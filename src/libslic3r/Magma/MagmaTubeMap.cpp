#include "MagmaTubeMap.hpp"

#include "../Layer.hpp"
#include "../Print.hpp"
#include "../Slicing.hpp"
#include "../ClipperUtils.hpp"
#include "../Polygon.hpp"
#include "../ExPolygon.hpp"
#include "../Surface.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
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
    if (layer_id < first_layer || layer_id > last_layer)
        return 0.0;
    int idx = layer_id - first_layer;
    if (idx < 0 || idx >= int(areas.size()))
        return 0.0;
    return areas[idx];
}

void CellPresence::mark_present(int layer_id, double area_val)
{
    if (first_layer == INT_MAX) {
        // First time: initialize
        first_layer = layer_id;
        last_layer  = layer_id;
        layers.assign(1, true);
        areas.assign(1, area_val);
        return;
    }

    // Expand range if needed
    if (layer_id < first_layer) {
        int extend = first_layer - layer_id;
        layers.insert(layers.begin(), extend, false);
        areas.insert(areas.begin(), extend, 0.0);
        first_layer = layer_id;
    }
    if (layer_id > last_layer) {
        int extend = layer_id - last_layer;
        layers.insert(layers.end(), extend, false);
        areas.insert(areas.end(), extend, 0.0);
        last_layer = layer_id;
    }

    int idx = layer_id - first_layer;
    layers[idx] = true;
    areas[idx]  = area_val;
}

// ============================================================================
// MagmaTubeMap — Statistics
// ============================================================================

int MagmaTubeMap::num_default_pairs() const
{
    int count = 0;
    for (const auto &p : m_pairs)
        if (!p.is_salvaged) ++count;
    return count;
}

int MagmaTubeMap::num_salvaged_pairs() const
{
    int count = 0;
    for (const auto &p : m_pairs)
        if (p.is_salvaged) ++count;
    return count;
}

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
    if (layer_id >= 0 && layer_id < int(m_print_z.size()))
        return m_print_z[layer_id];
    return 0.0;
}

double MagmaTubeMap::layer_height_at(int layer_id) const
{
    if (layer_id >= 0 && layer_id < int(m_layer_heights.size()))
        return m_layer_heights[layer_id];
    return double(m_layer_height);
}

double MagmaTubeMap::span_height_mm(int start_layer, int end_layer) const
{
    if (start_layer < 0 || end_layer < start_layer)
        return 0.0;
    int sz = int(m_print_z.size());
    if (start_layer >= sz || end_layer >= sz)
        return 0.0;
    // print_z is cumulative (top of layer), so span height =
    // top of end_layer minus bottom of start_layer.
    return m_print_z[end_layer] - (m_print_z[start_layer] - m_layer_heights[start_layer]);
}

int MagmaTubeMap::layer_at_height_from(int start_layer, double target_mm) const
{
    if (start_layer < 0 || start_layer >= int(m_print_z.size()))
        return m_num_layers - 1;
    // Target z = bottom of start_layer + target_mm
    double target_z = (m_print_z[start_layer] - m_layer_heights[start_layer]) + target_mm;
    // Binary search: find first layer with print_z >= target_z
    auto it = Slic3r::lower_bound_by_predicate(
        m_print_z.begin(), m_print_z.end(),
        [target_z](double pz) { return pz < target_z; });
    if (it == m_print_z.end())
        return m_num_layers - 1;
    return static_cast<int>(it - m_print_z.begin());
}

// ============================================================================
// MagmaTubeMap — Build
// ============================================================================

// Get effective interior width: use config value if > 0, otherwise auto-calculate
static float get_effective_interior_width(float config_value, float nozzle_diameter)
{
    return (config_value > 0) ? config_value
                              : static_cast<float>(calculate_auto_interior_width(nozzle_diameter));
}

std::unique_ptr<MagmaTubeMap> MagmaTubeMap::build(
    const std::vector<Layer*> &layers,
    const PrintRegionConfig &config,
    const PrintConfig &print_config,
    const PrintObjectConfig &obj_config,
    const SlicingParameters &slicing_params)
{
    auto t_start = std::chrono::high_resolution_clock::now();

    auto map = std::unique_ptr<MagmaTubeMap>(new MagmaTubeMap());

    // Extract config
    const float nozzle_diameter = static_cast<float>(print_config.nozzle_diameter.get_at(0));
    map->m_interior_width = get_effective_interior_width(
        config.magma_interior_width.value, nozzle_diameter);
    map->m_line_width = static_cast<float>(config.sparse_infill_line_width.get_abs_value(nozzle_diameter));
    if (map->m_line_width <= 0.f)
        map->m_line_width = nozzle_diameter;
    map->m_cell_spacing = cell_spacing_from_geometry(map->m_interior_width, map->m_line_width);
    // Nominal config layer height (fallback for missing layers in lookup tables)
    map->m_layer_height = static_cast<float>(obj_config.layer_height.value);
    map->m_dual_infill_enabled = config.dual_infill_enabled.value;

    // Build per-layer height/z tables for adaptive layer height support.
    // Must happen before WindowSpec and spiral params since they use m_min_layer_height.
    {
        int max_layer_id = 0;
        for (const Layer *l : layers)
            max_layer_id = std::max(max_layer_id, static_cast<int>(l->id()));
        map->m_num_layers = max_layer_id + 1;
        map->m_print_z.assign(map->m_num_layers, 0.0);
        map->m_layer_heights.assign(map->m_num_layers, 0.0);
        for (const Layer *l : layers) {
            int lid = static_cast<int>(l->id());
            map->m_print_z[lid]       = l->print_z;
            map->m_layer_heights[lid] = l->height;
        }
        // Use min_layer_height from SlicingParameters (already computed across
        // all extruders by OrcaSlicer). Clamp to nominal if something is off.
        map->m_min_layer_height = static_cast<float>(
            slicing_params.min_layer_height > 0
                ? std::min(slicing_params.min_layer_height, double(map->m_layer_height))
                : map->m_layer_height);
    }

    // Spiral params — use min_layer_height to cap helix angle at thinnest layer
    bool spiral_enabled = config.magma_spiral_interlock.value;
    map->m_spiral_params = compute_spiral_params(map->m_interior_width, map->m_line_width,
                                                  map->m_min_layer_height, spiral_enabled);

    // Window spec (handles window height and stagger; tube height may be overridden below).
    // Pass min_layer_height for conservative mm→layers conversions.
    map->m_window_spec = WindowSpec::from_config(
        config.magma_window_height.value,
        config.magma_tube_height_layers.value,
        static_cast<float>(config.magma_tube_height.value),
        map->m_interior_width,
        map->m_line_width,
        map->m_layer_height,
        config.magma_stagger_levels.value,
        map->m_min_layer_height
    );

    // Min tube height: window height (×2 for solid wall above+below) plus 2 layers of padding.
    // Use min_layer_height for the conversions so the mm limit is conservative.
    map->m_min_tube_height_mm = map->m_window_spec.window_height_layers * double(map->m_min_layer_height) * 2.0
                                + 2.0 * double(map->m_min_layer_height);
    map->m_min_tube_height_layers = map->m_window_spec.window_height_layers * 2 + 2;

    // -----------------------------------------------------------------------
    // Auto max tube height: coupled thermal + pressure model
    //
    // When the user hasn't specified an explicit tube height, we estimate
    // the maximum depth injected plastic can fill before either:
    //   (a) freezing in the channel (thermal limit), or
    //   (b) exceeding extruder force (pressure limit)
    //
    // These two constraints are coupled: deeper tubes require more pressure,
    // which means slower flow, which means more cooling time. Solving both
    // simultaneously eliminates volumetric speed as a variable.
    //
    // === Thermal constraint ===
    //
    //   Fill time:    t_fill  = 2 * depth * A_cell / V_dot
    //                 (U-tube path = 2x depth: down one cell, up the other)
    //   Freeze time:  t_freeze = D_h^2 / (4 * alpha)
    //                 (thermal diffusion across channel cross-section)
    //
    // === Pressure constraint (Hagen-Poiseuille, laminar flow) ===
    //
    //   Pressure drop: dP = 32 * mu * L * v / D_h^2
    //   With L = 2*depth (U-tube path) and v = V_dot / A_cell:
    //   V_dot_max = dP_max * A_cell * D_h^2 / (64 * mu * depth)
    //
    // === Coupled solution ===
    //
    //   Setting t_fill = fudge * t_freeze and V_dot = V_dot_max:
    //
    //     depth^2 = fudge * dP_max * D_h^4 / (512 * mu * alpha)
    //     depth   = D_h^2 * sqrt(fudge * dP_max / (512 * mu * alpha))
    //
    //   The volumetric speed drops out — it self-adjusts to balance both
    //   constraints. The implied injection speed is:
    //
    //     V_dot = dP_max * A_cell * D_h^2 / (64 * mu * depth)
    //
    // === Variables ===
    //
    //   A_cell = sqrt(3)/4 * iw^2     — triangular channel cross-section (mm^2)
    //   P_cell = 3 * iw               — triangular perimeter (mm)
    //   D_h    = 4 * A_cell / P_cell  — hydraulic diameter = 0.577 * iw (mm)
    //   alpha  ≈ 0.1 mm^2/s           — thermal diffusivity of FDM plastics
    //   mu     ≈ 50 Pa·s (= 50e-6 N·s/mm^2) — melt viscosity at injection
    //                                   temp, shear-thinned (~1000 s^-1)
    //   dP_max ≈ 15 MPa (= 15 N/mm^2) — typical direct-drive extruder max
    //                                   pressure (~40-60N on 1.75mm filament)
    //
    //   The fudge factor (magma_fill_depth_factor, default 3.0) absorbs:
    //     - Actual extruder force (varies by hardware)
    //     - Actual melt viscosity (varies by material and injection temp —
    //       higher injection temp significantly reduces viscosity, allowing
    //       deeper fills; increase fudge factor to account for this)
    //     - Wall insulation effects (surrounding plastic slows cooling)
    //     - Continuous flow bringing fresh heat into the channel
    //   Calibrate via test prints. Higher = taller tubes.
    //
    // === Optimality ===
    //
    //   The coupled solution is the global maximum depth. For any injection
    //   speed V: achievable depth = min(thermal(V), pressure(V)). Since
    //   thermal increases with V and pressure decreases with V, their
    //   intersection (the coupled solution) maximizes the min. No speed
    //   choice can produce deeper tubes.
    //
    // === Melt rate cap ===
    //
    //   The implied V_dot must not exceed filament_max_volumetric_speed,
    //   which represents the hotend's melt capacity. If it would, we reduce
    //   depth to the thermal-only limit at that capped speed (which is
    //   always the binding constraint for speeds below the coupled optimum).
    // -----------------------------------------------------------------------

    bool user_set_tube_height = (config.magma_tube_height_layers.value > 0
                                 || config.magma_tube_height.value > 0);

    if (!user_set_tube_height) {
        double iw = map->m_interior_width;
        double cell_area   = 0.433 * iw * iw;              // mm^2 (sqrt(3)/4 * iw^2)
        double cell_perim  = 3.0 * iw;                     // mm
        double hydraulic_d = 4.0 * cell_area / cell_perim;  // mm (= 0.577 * iw)
        double dh2 = hydraulic_d * hydraulic_d;

        constexpr double alpha   = 0.10;    // mm^2/s — thermal diffusivity
        constexpr double mu      = 50e-6;   // N·s/mm^2 (= 50 Pa·s) — melt viscosity
        constexpr double dp_max  = 15.0;    // N/mm^2 (= 15 MPa) — extruder max pressure

        double fudge = config.magma_fill_depth_factor.value;
        double max_depth_mm;
        double injection_vol_speed = config.magma_injection_speed.value;

        if (injection_vol_speed > 0) {
            // User specified injection speed — compute depth at that speed,
            // taking the tighter of thermal and pressure limits.
            //
            // Cap at hotend melt rate first: the actual G-code injection speed
            // is capped at max_vol, so tube depth must match achievable speed.
            // Without this cap, depth would be overly conservative (pressure
            // limit at a speed the extruder won't actually reach).
            double max_vol = print_config.filament_max_volumetric_speed.get_at(0);
            if (max_vol > 0 && injection_vol_speed > max_vol)
                injection_vol_speed = max_vol;

            // Thermal limit:  depth = fudge * V_dot * D_h^2 / (8 * A * alpha)
            // Pressure limit: depth = dP_max * A * D_h^2 / (64 * mu * V_dot)
            // For V < V_coupled: thermal is binding. For V > V_coupled: pressure is binding.
            double depth_thermal  = fudge * injection_vol_speed * dh2
                                    / (8.0 * cell_area * alpha);
            double depth_pressure = dp_max * cell_area * dh2
                                    / (64.0 * mu * injection_vol_speed);
            max_depth_mm = std::min(depth_thermal, depth_pressure);
        } else {
            // No explicit speed — use coupled pressure+thermal solution.
            // depth = D_h^2 * sqrt(fudge * dP_max / (512 * mu * alpha))
            max_depth_mm = dh2 * std::sqrt(fudge * dp_max / (512.0 * mu * alpha));

            // The coupled formula implies a specific injection speed:
            // V_dot = dP_max * A * D_h^2 / (64 * mu * depth)
            // Cap at hotend melt rate (filament_max_volumetric_speed) to avoid
            // exceeding what the hotend can melt. If capped, recompute using
            // thermal-only at the melt-limited speed.
            double implied_speed = dp_max * cell_area * dh2
                                   / (64.0 * mu * max_depth_mm);
            double max_vol = print_config.filament_max_volumetric_speed.get_at(0);
            if (max_vol > 0 && implied_speed > max_vol) {
                max_depth_mm = fudge * max_vol * dh2 / (8.0 * cell_area * alpha);
            }
        }

        map->m_max_tube_height_mm = std::max(max_depth_mm, map->m_min_tube_height_mm);
    } else {
        // User specified tube height — convert to mm (take the larger of layers and mm specs).
        double from_layers_mm = map->m_window_spec.tube_height_layers * double(map->m_layer_height);
        double from_mm = config.magma_tube_height.value;
        map->m_max_tube_height_mm = std::max({from_layers_mm, from_mm, map->m_min_tube_height_mm});
    }

    // Approximate layer count for logging and backward compat
    map->m_max_tube_height_layers = static_cast<int>(
        std::ceil(map->m_max_tube_height_mm / double(map->m_min_layer_height)));
    map->m_window_spec.tube_height_layers = map->m_max_tube_height_layers;

    // Build phases
    map->scan_layers(layers);
    map->detect_constrictions();
    map->assign_default_tubes();
    map->assign_salvage_tubes();
    map->compute_volumes(layers);

    auto t_end = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

    BOOST_LOG_TRIVIAL(info) << "MagmaTubeMap built in " << ms << "ms: "
        << map->num_cells() << " cells, "
        << map->num_default_pairs() << " default, "
        << map->num_salvaged_pairs() << " salvaged, "
        << map->num_solid_cells() << " solid | "
        << "tube_height=~" << map->m_max_tube_height_layers << "L ("
        << map->m_max_tube_height_mm << "mm"
        << (user_set_tube_height ? ", user" : ", auto") << ")";

    return map;
}

// ============================================================================
// MagmaTubeMap — scan_layers
// ============================================================================

void MagmaTubeMap::scan_layers(const std::vector<Layer*> &layers)
{
    // Pre-compute ideal inset triangle area (scaled^2) for interior cells.
    // Infill lines are centered on triangle edges, so the bead extends
    // line_width/2 inward from each edge. Polygon offset of an equilateral
    // triangle by -d reduces side by 2*d*sqrt(3). With d = line_width/2:
    //   inset_side = side - line_width * sqrt(3)
    const double half_line_width = m_line_width * 0.5;
    const double side = triangle_side_length(m_cell_spacing);
    const double inset_side = side - m_line_width * SQRT3;
    const double inset_area_mm2 = (inset_side > 0) ? 0.433 * inset_side * inset_side : 0.0;
    const double inset_area_scaled2 = inset_area_mm2 * 1e12;  // (1e6)^2

    // Boundary cells require ≥50% of ideal tube area to be unobstructed.
    // Below this, the tube cross-section is too constricted for plastic flow.
    const double min_area_scaled2 = inset_area_scaled2 * 0.50;

    // Interior inset: if cell center is this far inside the magma region,
    // the tube inscribed circle (diameter = interior_width) fits entirely.
    // This is the tube-flow-relevant check, not the triangle-corner check.
    const coord_t interior_inset = scale_(m_interior_width * 0.5);

    // FIXED reference lattice for cell IDENTITY (a,b,c coordinates).
    // Cell identity must be stable across layers for tube pairing to work.
    // Spiral offset is applied per-layer for POSITION checks below.
    TriangleLattice ref_lattice(m_cell_spacing, 0.0, 0.0);

    for (int i = 0; i < int(layers.size()); ++i) {
        const Layer *layer = layers[i];
        // Use Layer::id() — not the array index — so that pair_start_layer,
        // pair_end_layer, and CellPresence layer indices all match what
        // FillMagma and GCode will query with (Layer::id() includes raft offset).
        const int layer_id = static_cast<int>(layer->id());

        // Spiral-offset lattice for this layer's actual cell positions.
        // Cell identity comes from ref_lattice, but position checks use the
        // spiral-offset position (matching what FillMagma actually renders).
        Vec2d spiral_off = compute_spiral_offset(m_spiral_params, layer_id);
        TriangleLattice layer_lattice(m_cell_spacing, spiral_off.x(), spiral_off.y());

        // Collect surfaces where Magma infill will be generated
        ExPolygons magma_regions;
        for (const LayerRegion *layerm : layer->regions()) {
            const PrintRegionConfig &region_config = layerm->region().config();
            for (const Surface &surface : layerm->fill_surfaces.surfaces) {
                if (surface.is_zone_outer()) {
                    magma_regions.push_back(surface.expolygon);
                } else if (surface.surface_type == stInternal
                           && region_config.sparse_infill_pattern.value == ipMagmaTriangle) {
                    magma_regions.push_back(surface.expolygon);
                }
            }
        }

        magma_regions = union_ex(magma_regions);
        if (magma_regions.empty())
            continue;

        // Shrink region to identify fully-interior cells: tube inscribed circle
        // fits entirely within magma region → no per-cell clipping needed.
        ExPolygons interior_region = offset_ex(magma_regions, -interior_inset);

        // Enumerate cells from fixed reference lattice (stable identity)
        BoundingBox bbox = get_extents(magma_regions);
        // Expand bbox slightly to account for spiral offset moving cells
        bbox.offset(scale_(m_interior_width));
        std::vector<TriangleCell> cells = ref_lattice.enumerate_cells(bbox);

        for (const TriangleCell &cell : cells) {
            // Use spiral-offset position for containment checks
            Vec2d center_mm = layer_lattice.cell_center(cell);
            Point center_pt(scale_(center_mm.x()), scale_(center_mm.y()));

            // Fast path: center inside tube-clearance inset → fully unobstructed
            bool is_interior = false;
            for (const ExPolygon &ep : interior_region) {
                if (ep.contains(center_pt)) {
                    is_interior = true;
                    break;
                }
            }

            if (is_interior) {
                m_cells[cell].mark_present(layer_id, inset_area_scaled2);
                continue;
            }

            // Center must at least be inside the original region
            bool in_region = false;
            for (const ExPolygon &ep : magma_regions) {
                if (ep.contains(center_pt)) {
                    in_region = true;
                    break;
                }
            }

            if (!in_region)
                continue;

            // Boundary cell: compute actual tube area at spiral-offset position
            std::array<Vec2d, 3> corners = layer_lattice.cell_corners(cell);
            Polygon triangle;
            triangle.points.reserve(3);
            triangle.points.emplace_back(scale_(corners[0].x()), scale_(corners[0].y()));
            triangle.points.emplace_back(scale_(corners[1].x()), scale_(corners[1].y()));
            triangle.points.emplace_back(scale_(corners[2].x()), scale_(corners[2].y()));

            ExPolygons inset = offset_ex(triangle, -scale_(half_line_width));
            if (inset.empty())
                continue;

            ExPolygons clipped = intersection_ex(inset, magma_regions);
            if (clipped.empty())
                continue;

            double area = 0.0;
            for (const ExPolygon &ep : clipped)
                area += std::abs(ep.area());

            if (area < min_area_scaled2)
                continue;

            m_cells[cell].mark_present(layer_id, area);
        }
    }
}

// ============================================================================
// MagmaTubeMap — detect_constrictions
// ============================================================================

void MagmaTubeMap::detect_constrictions()
{
    // Two-tier constriction detection:
    // Tier 1 (fast): area ratio between consecutive layers
    // Tier 2 (expensive): polygon overlap check for flagged cells
    //
    // A constriction splits the cell's presence into separate spans,
    // preventing a tube from bridging across a near-discontinuity.

    int constrictions_found = 0;

    for (auto &[cell, presence] : m_cells) {
        if (presence.first_layer == presence.last_layer)
            continue;

        for (int i = presence.first_layer; i < presence.last_layer; ++i) {
            if (!presence.present(i) || !presence.present(i + 1))
                continue;

            double area_i   = presence.area(i);
            double area_ip1 = presence.area(i + 1);

            if (area_i <= 0.0 || area_ip1 <= 0.0)
                continue;

            // Tier 1: Area ratio heuristic (fast)
            double ratio = std::min(area_i, area_ip1) / std::max(area_i, area_ip1);
            if (ratio > 0.7)
                continue;  // Healthy overlap, skip expensive check

            // Tier 2: For cells that fail the ratio check, we use a stricter
            // area-based threshold. Full polygon overlap would require storing
            // per-layer magma regions (expensive). Instead, use a tighter ratio.
            if (ratio < 0.3) {
                // Severe constriction: split the cell at this boundary
                int idx = i + 1 - presence.first_layer;
                if (idx > 0 && idx < int(presence.layers.size())) {
                    presence.layers[idx] = false;
                    presence.areas[idx]  = 0.0;
                    ++constrictions_found;
                }
            }
        }
    }

    if (constrictions_found > 0)
        BOOST_LOG_TRIVIAL(info) << "MagmaTubeMap: " << constrictions_found << " constrictions found";
}

// ============================================================================
// MagmaTubeMap — assign_default_tubes
// ============================================================================

// Count shared present layers in a range
static int count_shared_layers(const CellPresence &a, const CellPresence &b,
                               int from, int to)
{
    int count = 0;
    for (int i = from; i <= to; ++i)
        if (a.present(i) && b.present(i))
            ++count;
    return count;
}

// Find contiguous shared spans (both cells present on every layer)
struct SharedSpan {
    int start;
    int end;  // inclusive
};

static std::vector<SharedSpan> find_shared_spans(const CellPresence &a, const CellPresence &b)
{
    std::vector<SharedSpan> spans;
    int range_start = std::max(a.first_layer, b.first_layer);
    int range_end   = std::min(a.last_layer,  b.last_layer);

    SharedSpan current{-1, -1};
    for (int i = range_start; i <= range_end; ++i) {
        if (a.present(i) && b.present(i)) {
            if (current.start < 0)
                current.start = i;
            current.end = i;
        } else {
            if (current.start >= 0) {
                spans.push_back(current);
                current = {-1, -1};
            }
        }
    }
    if (current.start >= 0)
        spans.push_back(current);

    return spans;
}

void MagmaTubeMap::assign_default_tubes()
{
    const double tube_h_mm = m_max_tube_height_mm;
    const int stagger_levels = m_window_spec.stagger_levels;

    // Only iterate up triangles — each up-down pair considered once.
    for (const auto &[cell, presence] : m_cells) {
        if (!cell.is_up())
            continue;

        // Try all 3 DN neighbors, pick the one with the longest total shared span.
        // This fixes poor coverage at sloped boundaries where the a-axis partner
        // only exists for a fraction of the cell's height but other neighbors
        // have much better overlap.
        auto neighbors = cell.neighbors();
        TriangleCell best_partner;
        double best_total_shared_mm = 0.0;
        const CellPresence *best_partner_presence = nullptr;

        for (const TriangleCell &neighbor : neighbors) {
            auto nit = m_cells.find(neighbor);
            if (nit == m_cells.end())
                continue;

            // Sum total shared height in mm across all contiguous spans
            std::vector<SharedSpan> nspans = find_shared_spans(presence, nit->second);
            double total_mm = 0.0;
            for (const SharedSpan &s : nspans)
                total_mm += span_height_mm(s.start, s.end);

            if (total_mm > best_total_shared_mm) {
                best_total_shared_mm = total_mm;
                best_partner = neighbor;
                best_partner_presence = &nit->second;
            }
        }

        if (!best_partner_presence || best_total_shared_mm < m_min_tube_height_mm)
            continue;

        const CellPresence &partner_presence = *best_partner_presence;
        TriangleCell partner = best_partner;

        // Stagger offset in mm: shifts tube grid so different cells have
        // tube tops/bottoms at different Z heights, avoiding weak planes.
        // Computed from UP cell coordinates — doesn't depend on partner choice.
        int stagger = cell.stagger_level(stagger_levels);
        double stagger_offset_mm = stagger * (tube_h_mm / stagger_levels);

        // Find contiguous shared spans with the best partner
        std::vector<SharedSpan> spans = find_shared_spans(presence, partner_presence);

        for (const SharedSpan &span : spans) {
            double span_h = span_height_mm(span.start, span.end);
            if (span_h < m_min_tube_height_mm)
                continue;

            // Find mm-based tube grid boundaries within this span.
            // The global grid has boundaries at z = stagger_offset_mm + i * tube_h_mm
            // (measured from z=0). We walk layers in the span, accumulating height,
            // and place a boundary when cumulative height crosses a grid line.
            double span_bottom_z = m_print_z[span.start] - m_layer_heights[span.start];

            // How far into the current tube_h_mm period is span_bottom_z?
            // Use fmod to find the offset, then first boundary is at that distance.
            double into_period = std::fmod(span_bottom_z - stagger_offset_mm, tube_h_mm);
            if (into_period < 0) into_period += tube_h_mm;
            double first_target_mm = tube_h_mm - into_period;
            if (first_target_mm < 1e-6) first_target_mm += tube_h_mm;  // on boundary exactly

            // Collect boundaries by walking layers and tracking cumulative height
            std::vector<int> boundaries;
            boundaries.push_back(span.start);

            double accum = 0.0;
            double target = first_target_mm;
            for (int i = span.start; i <= span.end; ++i) {
                accum += m_layer_heights[i];
                if (accum >= target - 1e-6 && i < span.end) {
                    // Place boundary AFTER this layer (start of next segment)
                    boundaries.push_back(i + 1);
                    target += tube_h_mm;
                }
            }
            boundaries.push_back(span.end + 1);  // one-past-end

            // Merge short first/last segments into their neighbors (mm-based check)
            while (boundaries.size() >= 3) {
                double first_seg_h = span_height_mm(boundaries[0], boundaries[1] - 1);
                if (first_seg_h >= m_min_tube_height_mm)
                    break;
                boundaries.erase(boundaries.begin() + 1);
            }
            while (boundaries.size() >= 3) {
                double last_seg_h = span_height_mm(
                    boundaries[boundaries.size() - 2], boundaries.back() - 1);
                if (last_seg_h >= m_min_tube_height_mm)
                    break;
                boundaries.erase(boundaries.end() - 2);
            }

            // Create tubes from segments
            for (size_t i = 0; i + 1 < boundaries.size(); ++i) {
                int seg_start = boundaries[i];
                int seg_end   = boundaries[i + 1] - 1;  // inclusive

                if (span_height_mm(seg_start, seg_end) < m_min_tube_height_mm)
                    continue;  // shouldn't happen after merging, but guard

                UTubePair pair;
                pair.cell_a = cell;
                pair.cell_b = partner;
                pair.pair_start_layer = seg_start;
                pair.pair_end_layer = seg_end;
                pair.is_salvaged = false;
                pair.volume_mm3 = 0.0;

                int pair_idx = static_cast<int>(m_pairs.size());
                m_pairs.push_back(pair);
                m_cell_pair_index[cell].push_back(pair_idx);
                m_cell_pair_index[partner].push_back(pair_idx);
            }
        }
    }
}

// ============================================================================
// MagmaTubeMap — assign_salvage_tubes
// ============================================================================

// Build a synthetic CellPresence that only marks layers NOT covered by
// existing tube pairs. Returns the uncovered layer ranges.
static CellPresence uncovered_presence(
    const CellPresence &presence,
    const std::vector<int> &pair_indices,
    const std::vector<UTubePair> &pairs)
{
    if (pair_indices.empty())
        return presence;

    CellPresence result;
    for (int layer = presence.first_layer; layer <= presence.last_layer; ++layer) {
        if (!presence.present(layer))
            continue;

        // Check if any existing pair covers this layer
        bool covered = false;
        for (int idx : pair_indices) {
            const UTubePair &p = pairs[idx];
            if (layer >= p.pair_start_layer && layer <= p.pair_end_layer) {
                covered = true;
                break;
            }
        }

        if (!covered)
            result.mark_present(layer, presence.area(layer));
    }

    return result;
}

void MagmaTubeMap::assign_salvage_tubes()
{
    // Collect cells with uncovered layer ranges, sorted by first_layer ascending.
    // A cell needs salvage if it has no entry in m_cell_pair_index (completely
    // unassigned) OR if it has existing pairs but some layers are uncovered.
    std::vector<TriangleCell> candidates;
    int completely_unassigned = 0;
    int partially_uncovered = 0;
    for (const auto &[cell, presence] : m_cells) {
        auto it = m_cell_pair_index.find(cell);
        if (it == m_cell_pair_index.end()) {
            // Completely unassigned
            candidates.push_back(cell);
            ++completely_unassigned;
        } else {
            // Check if there are uncovered layers
            CellPresence uc = uncovered_presence(presence, it->second, m_pairs);
            if (uc.first_layer <= uc.last_layer) {
                candidates.push_back(cell);
                ++partially_uncovered;
            }
        }
    }

    BOOST_LOG_TRIVIAL(info) << "Magma salvage: " << candidates.size() << " candidates ("
        << completely_unassigned << " unassigned, " << partially_uncovered << " partially uncovered) from "
        << m_cells.size() << " total cells, min_height=" << m_min_tube_height_mm << "mm"
        << ", max_height=" << m_max_tube_height_mm << "mm, total pairs before=" << m_pairs.size();

    std::sort(candidates.begin(), candidates.end(),
        [this](const TriangleCell &a, const TriangleCell &b) {
            return m_cells.at(a).first_layer < m_cells.at(b).first_layer;
        });

    int salvage_created = 0;
    int skip_already_covered = 0;
    int skip_no_neighbors = 0;
    int skip_too_short = 0;
    int fail_up = 0, fail_dn = 0;
    int max_uncov_wasted = 0;  // largest uncovered span that failed

    for (const TriangleCell &cell : candidates) {
        const CellPresence &presence = m_cells.at(cell);

        // Get this cell's uncovered layers
        auto cell_idx_it = m_cell_pair_index.find(cell);
        std::vector<int> cell_pairs;
        if (cell_idx_it != m_cell_pair_index.end())
            cell_pairs = cell_idx_it->second;

        CellPresence cell_uncovered = uncovered_presence(presence, cell_pairs, m_pairs);
        if (cell_uncovered.first_layer > cell_uncovered.last_layer) {
            ++skip_already_covered;
            continue;  // Fully covered (may have been paired by earlier iteration)
        }

        int uncovered_height = cell_uncovered.last_layer - cell_uncovered.first_layer + 1;

        // Try all 3 neighbors, prefer longest shared uncovered height (in mm)
        TriangleCell best_partner;
        double best_height_mm = 0.0;
        SharedSpan best_span{-1, -1};

        int neighbors_missing = 0;
        int neighbors_covered = 0;
        int neighbors_too_short = 0;

        for (const TriangleCell &neighbor : cell.neighbors()) {
            // Skip if not present in cell map
            auto nit = m_cells.find(neighbor);
            if (nit == m_cells.end()) {
                ++neighbors_missing;
                continue;
            }
            const CellPresence &neighbor_presence = nit->second;

            // Get neighbor's uncovered layers
            auto nbr_idx_it = m_cell_pair_index.find(neighbor);
            std::vector<int> nbr_pairs;
            if (nbr_idx_it != m_cell_pair_index.end())
                nbr_pairs = nbr_idx_it->second;

            CellPresence nbr_uncovered = uncovered_presence(neighbor_presence, nbr_pairs, m_pairs);
            if (nbr_uncovered.first_layer > nbr_uncovered.last_layer) {
                ++neighbors_covered;
                continue;  // Neighbor fully covered
            }

            // Find shared spans within uncovered ranges
            std::vector<SharedSpan> spans = find_shared_spans(cell_uncovered, nbr_uncovered);

            bool any_tall_enough = false;
            for (const SharedSpan &span : spans) {
                double h_mm = span_height_mm(span.start, span.end);
                if (h_mm >= m_min_tube_height_mm && h_mm > best_height_mm) {
                    best_partner = neighbor;
                    best_height_mm = h_mm;
                    best_span = span;
                    any_tall_enough = true;
                }
            }
            if (!any_tall_enough)
                ++neighbors_too_short;
        }

        if (best_height_mm > 0) {
            // Cap to max tube height (mm-based)
            int seg_end = std::min(
                layer_at_height_from(best_span.start, m_max_tube_height_mm),
                best_span.end);

            UTubePair pair;
            pair.cell_a = cell;
            pair.cell_b = best_partner;
            pair.pair_start_layer = best_span.start;
            pair.pair_end_layer = seg_end;
            pair.is_salvaged = true;
            pair.volume_mm3 = 0.0;

            int pair_idx = static_cast<int>(m_pairs.size());
            m_pairs.push_back(pair);
            m_cell_pair_index[cell].push_back(pair_idx);
            m_cell_pair_index[best_partner].push_back(pair_idx);
            ++salvage_created;
        } else {
            // No partner found — mark as solid fill (empty vector)
            if (m_cell_pair_index.find(cell) == m_cell_pair_index.end())
                m_cell_pair_index[cell];  // default-inserts empty vector

            if (neighbors_missing == 3)
                ++skip_no_neighbors;
            else
                ++skip_too_short;

            if (cell.is_up()) ++fail_up; else ++fail_dn;
            max_uncov_wasted = std::max(max_uncov_wasted, uncovered_height);
        }
    }

    BOOST_LOG_TRIVIAL(info) << "Magma salvage results: " << salvage_created << " created, "
        << skip_already_covered << " covered, " << skip_no_neighbors << " no_nbrs, "
        << skip_too_short << " too_short (min=" << m_min_tube_height_mm << "mm) | fails: "
        << fail_up << " UP " << fail_dn << " DN, max_wasted=" << max_uncov_wasted << "L"
        << " | total pairs now=" << m_pairs.size();
}

// ============================================================================
// MagmaTubeMap — compute_volumes
// ============================================================================

void MagmaTubeMap::compute_volumes(const std::vector<Layer*> &layers)
{
    // Build layer_id → height map for per-layer volume calculation.
    // Handles initial layer height, adaptive/dynamic layer heights, and raft offsets.
    std::unordered_map<int, double> layer_heights;
    for (const Layer *layer : layers)
        layer_heights[static_cast<int>(layer->id())] = layer->height;

    double edge_len = triangle_side_length(m_cell_spacing);

    for (UTubePair &pair : m_pairs) {
        double tube_volume_scaled2_mm = 0.0;  // accumulated (area_scaled2 * height_mm)
        double window_height_mm = 0.0;

        for (int layer_id = pair.pair_start_layer; layer_id <= pair.pair_end_layer; ++layer_id) {
            auto it_a = m_cells.find(pair.cell_a);
            auto it_b = m_cells.find(pair.cell_b);

            double area_a = (it_a != m_cells.end()) ? it_a->second.area(layer_id) : 0.0;
            double area_b = (it_b != m_cells.end()) ? it_b->second.area(layer_id) : 0.0;

            auto h_it = layer_heights.find(layer_id);
            double lh = (h_it != layer_heights.end()) ? h_it->second : double(m_layer_height);

            tube_volume_scaled2_mm += (area_a + area_b) * lh;

            // Accumulate window height for window layers
            if (layer_id >= pair.pair_start_layer &&
                layer_id < pair.pair_start_layer + m_window_spec.window_height_layers)
                window_height_mm += lh;
        }

        // Convert: area in scaled^2 → mm^2 via SCALING_FACTOR^2 (1e-12)
        double tube_volume = tube_volume_scaled2_mm * SCALING_FACTOR * SCALING_FACTOR;

        // Window gap volume: shared edge opening between paired cells.
        double window_volume = edge_len * m_line_width * window_height_mm;

        pair.volume_mm3 = tube_volume + window_volume;
    }
}

// ============================================================================
// MagmaTubeMap — Query Interface
// ============================================================================

bool MagmaTubeMap::is_window_open(const TriangleCell &cell, int layer_id) const
{
    auto it = m_cell_pair_index.find(cell);
    if (it == m_cell_pair_index.end() || it->second.empty())
        return false;

    // Windows are at the bottom of each tube segment — the first
    // window_height_layers layers. This is where the U-tube connection
    // gap allows plastic to flow from one cell to the other.
    // Injection enters at pair_end_layer (top) and flows down.
    const int wh = m_window_spec.window_height_layers;

    for (int pair_idx : it->second) {
        const UTubePair &pair = m_pairs[pair_idx];
        if (layer_id >= pair.pair_start_layer &&
            layer_id < pair.pair_start_layer + wh)
            return true;
    }

    return false;
}

bool MagmaTubeMap::is_paired(const TriangleCell &cell) const
{
    auto it = m_cell_pair_index.find(cell);
    return (it != m_cell_pair_index.end() && !it->second.empty());
}

// Merge overlapping/adjacent intervals in a sorted vector.
static void merge_intervals(std::vector<std::pair<double, double>> &intervals)
{
    if (intervals.size() <= 1)
        return;
    std::sort(intervals.begin(), intervals.end());
    std::vector<std::pair<double, double>> merged;
    merged.push_back(intervals[0]);
    for (size_t i = 1; i < intervals.size(); ++i) {
        if (intervals[i].first <= merged.back().second + 0.01)
            merged.back().second = std::max(merged.back().second, intervals[i].second);
        else
            merged.push_back(intervals[i]);
    }
    intervals = std::move(merged);
}

WindowGaps MagmaTubeMap::window_gaps(int layer_id) const
{
    Vec2d offset = compute_spiral_offset(m_spiral_params, layer_id);
    TriangleLattice lattice(m_cell_spacing, offset.x(), offset.y());

    WindowGaps result;
    result.cell_spacing = m_cell_spacing;
    result.offset_x = offset.x();
    result.offset_y = offset.y();

    // For each pair with an open window on this layer, determine the shared
    // edge between cell_a and cell_b and add a gap on the correct line family.
    for (const auto &pair : m_pairs) {
        if (!is_window_open(pair.cell_a, layer_id))
            continue;

        SharedEdge edge = shared_edge(pair.cell_a, pair.cell_b);

        switch (edge) {
        case SharedEdge::Horizontal: {
            // Shared edge at row = max(cell_a.b, cell_b.b) for UP cell's base,
            // or min(cell_a.b, cell_b.b) + 1 for DOWN→UP boundary.
            // The shared vertices are (a, b) and (a+1, b) where b is the
            // UP cell's b coordinate.
            const TriangleCell &up = pair.cell_a.is_up() ? pair.cell_a : pair.cell_b;
            int row = up.b;
            Vec2d v0 = lattice.to_world(up.a, row);
            Vec2d v1 = lattice.to_world(up.a + 1, row);
            double x_lo = std::min(v0.x(), v1.x());
            double x_hi = std::max(v0.x(), v1.x());
            result.horiz[row].push_back({x_lo, x_hi});
            break;
        }
        case SharedEdge::Col60: {
            // Shared edge on column a (constant lattice x).
            // The shared vertices are (a, b) and (a, b+1) where a = max of the
            // two cells' a coordinates.
            int col = std::max(pair.cell_a.a, pair.cell_b.a);
            int row = std::min(pair.cell_a.b, pair.cell_b.b);
            Vec2d v0 = lattice.to_world(col, row);
            Vec2d v1 = lattice.to_world(col, row + 1);
            double y_lo = std::min(v0.y(), v1.y());
            double y_hi = std::max(v0.y(), v1.y());
            result.col60[col].push_back({y_lo, y_hi});
            break;
        }
        case SharedEdge::Diag120: {
            // Shared edge on diagonal s = a + b + 1 of the UP cell.
            // Shared vertices: UP's (a+1, b) and (a, b+1).
            const TriangleCell &up = pair.cell_a.is_up() ? pair.cell_a : pair.cell_b;
            int diag = up.a + up.b + 1;
            Vec2d v0 = lattice.to_world(up.a + 1, up.b);
            Vec2d v1 = lattice.to_world(up.a, up.b + 1);
            double y_lo = std::min(v0.y(), v1.y());
            double y_hi = std::max(v0.y(), v1.y());
            result.diag120[diag].push_back({y_lo, y_hi});
            break;
        }
        }
    }

    // Sort and merge intervals per line index
    for (auto &[key, intervals] : result.horiz)
        merge_intervals(intervals);
    for (auto &[key, intervals] : result.col60)
        merge_intervals(intervals);
    for (auto &[key, intervals] : result.diag120)
        merge_intervals(intervals);

    return result;
}

} // namespace magma
} // namespace Slic3r
