#include "MagmaTubeMap.hpp"

#include "../Layer.hpp"
#include "../Print.hpp"
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
    const PrintObjectConfig &obj_config)
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
    map->m_layer_height = layers.empty() ? 0.2f : static_cast<float>(layers.front()->height);
    map->m_num_layers = static_cast<int>(layers.size());
    map->m_dual_infill_enabled = config.dual_infill_enabled.value;

    // Spiral params
    bool spiral_enabled = config.magma_spiral_interlock.value;
    map->m_spiral_params = compute_spiral_params(map->m_interior_width, map->m_line_width, map->m_layer_height, spiral_enabled);

    // Window spec (handles window height and stagger; tube height may be overridden below)
    map->m_window_spec = WindowSpec::from_config(
        config.magma_window_height.value,
        config.magma_tube_height_layers.value,
        static_cast<float>(config.magma_tube_height.value),
        map->m_interior_width,
        map->m_line_width,
        map->m_layer_height,
        config.magma_stagger_levels.value
    );

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

        int auto_max_layers = static_cast<int>(std::ceil(max_depth_mm / map->m_layer_height));
        auto_max_layers = std::max(auto_max_layers, map->m_min_tube_height_layers);

        map->m_window_spec.tube_height_layers = auto_max_layers;
    }

    map->m_max_tube_height_layers = map->m_window_spec.tube_height_layers;

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
        << "tube_height=" << map->m_max_tube_height_layers << "L ("
        << map->m_max_tube_height_layers * (double)map->m_layer_height << "mm"
        << (user_set_tube_height ? ", user" : ", auto") << ")";

    return map;
}

// ============================================================================
// MagmaTubeMap — scan_layers
// ============================================================================

void MagmaTubeMap::scan_layers(const std::vector<Layer*> &layers)
{
    // Pre-compute ideal inset triangle area (scaled^2) for interior cells.
    // Interior cells (center fully inside magma region) use this constant;
    // only boundary cells need expensive polygon clipping.
    const double half_line_width = m_line_width * 0.5;
    const double side = triangle_side_length(m_cell_spacing);
    const double ideal_area_mm2 = 0.433 * side * side;  // sqrt(3)/4 * side^2
    // Inset reduces effective side by ~line_width * 2/sqrt(3)
    const double inset_side = side - m_line_width * 2.0 * INV_SQRT3;
    const double inset_area_mm2 = (inset_side > 0) ? 0.433 * inset_side * inset_side : 0.0;
    const double inset_area_scaled2 = inset_area_mm2 * 1e12;  // (1e6)^2

    const double min_area_mm2 = ideal_area_mm2 * 0.10;
    const double min_area_scaled2 = min_area_mm2 * 1e12;

    // Inset distance for shrinking magma region to test "fully interior" cells.
    // A cell center inside the shrunk region means ALL 3 corners are inside the original.
    const coord_t interior_inset = scale_(m_cell_spacing * 0.6);  // conservative: > incircle radius

    // IMPORTANT: Use a FIXED reference lattice (no spiral offset) for cell identity.
    // The spiral offset shifts the lattice origin per layer, which would produce
    // different (a,b,c) coordinates for the same physical cell on each layer.
    // Cell identity must be stable across layers for tube pairing to work.
    // The spiral offset is only applied during infill rendering (FillMagma.cpp).
    TriangleLattice ref_lattice(m_cell_spacing, 0.0, 0.0);

    for (int i = 0; i < int(layers.size()); ++i) {
        const Layer *layer = layers[i];
        // Use Layer::id() — not the array index — so that pair_start_layer,
        // pair_end_layer, and CellPresence layer indices all match what
        // FillMagma and GCode will query with (Layer::id() includes raft offset).
        const int layer_id = static_cast<int>(layer->id());

        // Collect surfaces where Magma infill will be generated
        ExPolygons magma_regions;
        for (const LayerRegion *layerm : layer->regions()) {
            const PrintRegionConfig &region_config = layerm->region().config();
            for (const Surface &surface : layerm->fill_surfaces.surfaces) {
                if (surface.is_zone_outer()) {
                    magma_regions.push_back(surface.expolygon);
                } else if (surface.surface_type == stInternal
                           && !m_dual_infill_enabled
                           && region_config.sparse_infill_pattern.value == ipMagmaTriangle) {
                    magma_regions.push_back(surface.expolygon);
                }
            }
        }

        magma_regions = union_ex(magma_regions);
        if (magma_regions.empty())
            continue;

        // Shrink region once to identify fully-interior cells (no per-cell clipping needed)
        ExPolygons interior_region = offset_ex(magma_regions, -interior_inset);

        // Enumerate cells from fixed reference lattice covering zone bbox
        BoundingBox bbox = get_extents(magma_regions);
        std::vector<TriangleCell> cells = ref_lattice.enumerate_cells(bbox);

        for (const TriangleCell &cell : cells) {
            Vec2d center_mm = ref_lattice.cell_center(cell);
            Point center_pt(scale_(center_mm.x()), scale_(center_mm.y()));

            // Fast path: cell center inside shrunk region → fully interior
            bool is_interior = false;
            for (const ExPolygon &ep : interior_region) {
                if (ep.contains(center_pt)) {
                    is_interior = true;
                    break;
                }
            }

            if (is_interior) {
                // Interior cell: use pre-computed ideal inset area
                m_cells[cell].mark_present(layer_id, inset_area_scaled2);
                continue;
            }

            // Check if center is at least inside the original region
            bool in_region = false;
            for (const ExPolygon &ep : magma_regions) {
                if (ep.contains(center_pt)) {
                    in_region = true;
                    break;
                }
            }

            if (!in_region)
                continue;  // Cell center outside magma region entirely

            // Boundary cell: do the expensive clip to get actual area
            std::array<Vec2d, 3> corners = ref_lattice.cell_corners(cell);
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
    const int tube_h = m_max_tube_height_layers;
    const int stagger_levels = m_window_spec.stagger_levels;

    // Only iterate up triangles — each up-down pair considered once.
    for (const auto &[cell, presence] : m_cells) {
        if (!cell.is_up())
            continue;

        TriangleCell partner = cell.get_paired_cell();  // a-axis down neighbor

        auto it = m_cells.find(partner);
        if (it == m_cells.end())
            continue;

        const CellPresence &partner_presence = it->second;

        // Phase offset: shifts tube boundaries so different cells have
        // tube tops/bottoms at different layers, avoiding weak planes.
        int stagger = cell.stagger_level(stagger_levels);
        int phase_offset = stagger * (tube_h / stagger_levels);

        // Find contiguous shared spans
        std::vector<SharedSpan> spans = find_shared_spans(presence, partner_presence);

        for (const SharedSpan &span : spans) {
            int shared_height = span.end - span.start + 1;

            if (shared_height < m_min_tube_height_layers)
                continue;

            // Find the phase-shifted tube grid boundaries within this span.
            // Grid boundaries occur at layers: phase_offset + i * tube_h (for integer i).
            // We find the first boundary at or after span.start, then enumerate segments.

            // First grid boundary >= span.start
            int first_boundary;
            if (phase_offset >= span.start) {
                first_boundary = phase_offset;
            } else {
                // Find smallest (phase_offset + i*tube_h) >= span.start
                int steps = (span.start - phase_offset + tube_h - 1) / tube_h;
                first_boundary = phase_offset + steps * tube_h;
            }

            // Collect segment boundaries: [span.start, boundary1, boundary2, ..., span.end+1)
            std::vector<int> boundaries;
            boundaries.push_back(span.start);
            for (int b = first_boundary; b <= span.end; b += tube_h) {
                if (b > span.start && b <= span.end)
                    boundaries.push_back(b);
            }
            boundaries.push_back(span.end + 1);  // one-past-end

            // Merge short first/last segments into their neighbors
            // A segment is "short" if it's less than m_min_tube_height_layers
            while (boundaries.size() >= 3) {
                int first_seg = boundaries[1] - boundaries[0];
                if (first_seg >= m_min_tube_height_layers)
                    break;
                // Merge first segment into second by removing the first interior boundary
                boundaries.erase(boundaries.begin() + 1);
            }
            while (boundaries.size() >= 3) {
                int last_seg = boundaries.back() - boundaries[boundaries.size() - 2];
                if (last_seg >= m_min_tube_height_layers)
                    break;
                // Merge last segment into second-to-last
                boundaries.erase(boundaries.end() - 2);
            }

            // Create tubes from segments
            for (size_t i = 0; i + 1 < boundaries.size(); ++i) {
                int seg_start = boundaries[i];
                int seg_end   = boundaries[i + 1] - 1;  // inclusive
                int seg_height = seg_end - seg_start + 1;

                if (seg_height < m_min_tube_height_layers)
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
    for (const auto &[cell, presence] : m_cells) {
        auto it = m_cell_pair_index.find(cell);
        if (it == m_cell_pair_index.end()) {
            // Completely unassigned
            candidates.push_back(cell);
        } else {
            // Check if there are uncovered layers
            CellPresence uc = uncovered_presence(presence, it->second, m_pairs);
            if (uc.first_layer <= uc.last_layer)
                candidates.push_back(cell);
        }
    }

    std::sort(candidates.begin(), candidates.end(),
        [this](const TriangleCell &a, const TriangleCell &b) {
            return m_cells.at(a).first_layer < m_cells.at(b).first_layer;
        });

    for (const TriangleCell &cell : candidates) {
        const CellPresence &presence = m_cells.at(cell);

        // Get this cell's uncovered layers
        auto cell_idx_it = m_cell_pair_index.find(cell);
        std::vector<int> cell_pairs;
        if (cell_idx_it != m_cell_pair_index.end())
            cell_pairs = cell_idx_it->second;

        CellPresence cell_uncovered = uncovered_presence(presence, cell_pairs, m_pairs);
        if (cell_uncovered.first_layer > cell_uncovered.last_layer)
            continue;  // Fully covered (may have been paired by earlier iteration)

        // Try all 3 neighbors, prefer longest shared uncovered height
        TriangleCell best_partner;
        int best_height = 0;
        SharedSpan best_span{-1, -1};

        for (const TriangleCell &neighbor : cell.neighbors()) {
            // Skip if not present in cell map
            auto nit = m_cells.find(neighbor);
            if (nit == m_cells.end())
                continue;

            const CellPresence &neighbor_presence = nit->second;

            // Get neighbor's uncovered layers
            auto nbr_idx_it = m_cell_pair_index.find(neighbor);
            std::vector<int> nbr_pairs;
            if (nbr_idx_it != m_cell_pair_index.end())
                nbr_pairs = nbr_idx_it->second;

            CellPresence nbr_uncovered = uncovered_presence(neighbor_presence, nbr_pairs, m_pairs);
            if (nbr_uncovered.first_layer > nbr_uncovered.last_layer)
                continue;  // Neighbor fully covered

            // Find shared spans within uncovered ranges
            std::vector<SharedSpan> spans = find_shared_spans(cell_uncovered, nbr_uncovered);

            for (const SharedSpan &span : spans) {
                int height = span.end - span.start + 1;
                if (height >= m_min_tube_height_layers && height > best_height) {
                    best_partner = neighbor;
                    best_height = height;
                    best_span = span;
                }
            }
        }

        if (best_height > 0) {
            // Cap to max tube height
            int seg_end = std::min(best_span.start + m_max_tube_height_layers - 1, best_span.end);

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
        } else {
            // No partner found — mark as solid fill (empty vector)
            if (m_cell_pair_index.find(cell) == m_cell_pair_index.end())
                m_cell_pair_index[cell];  // default-inserts empty vector
        }
    }
}

// ============================================================================
// MagmaTubeMap — compute_volumes
// ============================================================================

void MagmaTubeMap::compute_volumes(const std::vector<Layer*> &layers)
{
    for (UTubePair &pair : m_pairs) {
        double total_area_scaled2 = 0.0;

        for (int layer_id = pair.pair_start_layer; layer_id <= pair.pair_end_layer; ++layer_id) {
            auto it_a = m_cells.find(pair.cell_a);
            auto it_b = m_cells.find(pair.cell_b);

            double area_a = (it_a != m_cells.end()) ? it_a->second.area(layer_id) : 0.0;
            double area_b = (it_b != m_cells.end()) ? it_b->second.area(layer_id) : 0.0;

            total_area_scaled2 += area_a + area_b;
        }

        // Convert from scaled^2 * layers to mm^3
        // area in scaled^2 → mm^2: divide by 1e12 (SCALING_FACTOR^2)
        // volume = area_mm2 * layer_height
        double area_mm2 = total_area_scaled2 * SCALING_FACTOR * SCALING_FACTOR;
        pair.volume_mm3 = area_mm2 * m_layer_height;
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
