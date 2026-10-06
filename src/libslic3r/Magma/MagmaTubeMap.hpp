#ifndef slic3r_Magma_MagmaTubeMap_hpp_
#define slic3r_Magma_MagmaTubeMap_hpp_

#include "MagmaTriangleCell.hpp"
#include "MagmaSpiralOffset.hpp"
#include "../PrintConfig.hpp"
#include "../ExPolygon.hpp"
#include "../BoundingBox.hpp"

#include <functional>
#include <map>
#include <memory>
#include <unordered_map>
#include <optional>
#include <vector>
#include <climits>

namespace Slic3r {

class Layer;
struct SlicingParameters;

namespace magma {

// The tube interior width is configured directly (magma_interior_width), floored at 0.1mm.
inline double effective_interior_width(double interior_width)
{
    return std::max(0.1, interior_width);
}

// The pattern that prints the lattice: the outer zone pattern with dual infill, else the sparse
// infill pattern. It may be a non-Magma pattern; callers guard with is_magma_pattern().
inline InfillPattern magma_effective_pattern(const PrintRegionConfig &config) {
    return config.dual_infill_enabled.value ? config.dual_infill_outer_pattern.value
                                            : config.sparse_infill_pattern.value;
}

// Per-cell layer presence and interior area
struct CellPresence {
    int first_layer = INT_MAX;   // first layer where cell exists in magma region
    int last_layer  = -1;        // last layer where cell exists in magma region
    std::vector<bool>   layers;  // indexed [layer_id - first_layer]
    std::vector<double> areas;   // parallel, interior area in scaled^2 units
    std::vector<double> distances; // parallel, distance to nearest boundary (mm, unscaled)
    std::vector<double> opening_radii; // parallel, max injection-point→opening boundary (mm)
    std::vector<Vec2d>  injection_pts; // parallel, injection point (mm): clipped-opening
                                       // centroid; == cell centre for unclipped cells

    bool   present(int layer_id) const;
    // These return data only where the cell is present. A cell can be absent on layers inside
    // [first_layer, last_layer] (a hole passing through it); those rows hold placeholder zeros,
    // so injection_point() returns nullopt there rather than the world origin.
    double area(int layer_id) const;
    double distance(int layer_id) const;
    double opening_radius(int layer_id) const;
    std::optional<Vec2d> injection_point(int layer_id) const;
    void   mark_present(int layer_id, double area_val, double dist_val,
                        double opening_r, const Vec2d &inj_pt);
};

// A U-tube pair assignment — for tri-hex, a manifold (one hub + N equal-length vent
// legs). cell_a is the injection HUB, cell_b the PRIMARY vent leg (the 1-leg U-tube
// for triangle/square). extra_vents holds any additional legs added by the post-solver
// extra-vent sweep (empty for triangle/square). All legs span the same [start,cap] with
// windows pinned to the tube bottom; one injection fills the hub + every leg.
struct UTubePair {
    CellId cell_a;          // injection side (HUB for tri-hex)
    CellId cell_b;          // primary vent leg
    std::vector<CellId> extra_vents;  // tri-hex manifold: additional vent legs (sweep-assigned)
    int  pair_start_layer;        // first layer of this tube segment
    int  pair_end_layer;          // last layer of this tube segment (inclusive)
    double window_end_z;          // Z coordinate (mm) where the window ends — pure mm, no rounding
    double volume_mm3;            // injection volume (hub + all legs combined)

    // injection_center starts as the scan's clipped-cell centroid; measure_volumes() refines it
    // to the deposited cap-cavity centroid and sets measured_cap_opening_dia from the same
    // cavity, so center and seal opening always agree.
    Vec2d  injection_center;      // XY center for injection
    int    window_center_layer;   // center layer of window gap
    double measured_cap_opening_dia = 0.0;  // real cap opening (0 → cap_opening_diameter falls
                                            // back to the cell_inset model)
};

// Per-layer data: Z heights and pre-built lattice with spiral offset.
struct LayerData {
    double print_z;   // cumulative Z (top of layer)
    double height;    // individual layer height
    // Lattice with this layer's spiral offset. shared_ptr because LayerData is copied and
    // MagmaLattice is abstract.
    std::shared_ptr<MagmaLattice> lattice;

    double bottom_z() const { return print_z - height; }
};

class MagmaTubeMap {
public:
    // Build the tube map. Called from PrintObject::infill() before TBB fill loop.
    using ProgressFn      = std::function<void(int, int)>; // (current, total)
    using ThrowIfCanceled = std::function<void()>;

    static std::unique_ptr<MagmaTubeMap> build(
        const std::vector<Layer*> &layers,
        const PrintRegionConfig &config,
        const PrintConfig &print_config,
        const PrintObjectConfig &obj_config,
        const SlicingParameters &slicing_params,
        ProgressFn progress_fn = nullptr,
        ThrowIfCanceled throw_if_canceled = nullptr);

    // === Query interface (thread-safe, const) ===

    // Is this pair's window (wall gap) open on this layer? Z-only; the per-shape toolpaths
    // compute the cut geometry themselves.
    bool window_open_at(const UTubePair &pair, int layer_id) const;

    // Is this cell assigned to a tube (true) or should be solid fill (false)?
    bool is_paired(const CellId &cell) const;

    // Accessors for shared params (so FillMagma uses identical values)
    const SpiralParams& spiral_params() const { return m_spiral_params; }
    double cell_spacing() const { return m_cell_spacing; }
    float  interior_width() const { return m_interior_width; }

    // Nominal opening the nozzle must cover to seal: circumscribed diameter of the unclipped
    // inset cell.
    double tube_opening_diameter() const;
    // Opening at a pair's cap layer, for the per-tube seal depth: the measured cap cavity if
    // measure_volumes() has run, else the scan's clipped opening, else tube_opening_diameter().
    double cap_opening_diameter(const UTubePair &pair) const;
    // Preview-only bore (mm) of a cell's injected column: the ideal bore scaled by
    // sqrt(area(layer) / max_area), so the column narrows where the part clips the cell.
    double cell_bore_at(const CellId &cell, int layer) const;

    // Layers below first_layer_id() (raft layers) carry no data or lattice.
    bool   valid_layer(int layer_id) const {
        return layer_id >= m_first_layer_id && layer_id < int(m_layer_data.size());
    }
    int    first_layer_id() const { return m_first_layer_id; }
    // The pattern this map was built for. One map serves the whole object, so Fill.cpp rejects
    // regions that select a different Magma pattern.
    InfillPattern pattern() const { return m_pattern; }
    // Lattice with this layer's spiral offset applied.
    const MagmaLattice& lattice_at(int layer_id) const { return *m_layer_data[layer_id].lattice; }

    // First layer at or above the Z midpoint of the pair's window gap.
    int window_center_layer(const UTubePair& pair) const;

    // For injection G-code and visualization
    const std::vector<UTubePair>& u_tube_pairs() const { return m_pairs; }

    // Measure injected volumes from the deposited toolpaths; call after PrintObject::infill().
    // Sets pair.volume_mm3 and rebuilds the injection cap-layer list.
    void measure_volumes(const std::vector<Layer*> &layers);

    float layer_height() const { return m_layer_height; }
    double window_height_mm() const { return m_window_spec.window_height_mm; }

    // Per-layer data accessors (adaptive layer height support)
    double print_z(int layer_id) const;
    double layer_height_at(int layer_id) const;

    // Slicing warnings for the user, one per independent condition found during build.
    const std::vector<std::string>& warning_messages() const { return m_warning_messages; }

    // Layer IDs that have injection caps (pair_end_layer of each UTubePair).
    // Used by ToolOrdering to register injection filament on the correct layers.
    const std::vector<int>& injection_layer_ids() const { return m_injection_layer_ids; }

    // 0-based extruder that performs the injections, resolved once at build time. The cap
    // clearances here are planned for its nozzle, so the emitter and ToolOrdering read it back
    // instead of resolving it again.
    int injection_extruder() const { return m_injection_extruder; }

    // Statistics for debug logging
    int num_cells() const { return static_cast<int>(m_cells.size()); }
    int num_pairs() const { return static_cast<int>(m_pairs.size()); }
    int num_solid_cells() const;

private:
    MagmaTubeMap() = default;

    // Build phases
    void scan_layers(const std::vector<Layer*> &layers);
    void assign_tubes(ProgressFn progress_fn, ThrowIfCanceled throw_if_canceled);
    // Tri-hex only: gives vents extra legs on bordering hub-tubes after the solver.
    void assign_extra_vents();
    void precompute_window_end_z();
    void precompute_injection_data();
    // Lazily built zero-offset lattice for cell identity and topology, which do not depend on
    // the spiral offset.
    const MagmaLattice& topology_lattice() const;

    // Adaptive layer height helpers
    double span_height_mm(int start_layer, int end_layer) const;

    // Data
    std::unordered_map<CellId, CellPresence, CellIdHash> m_cells;
    std::vector<UTubePair> m_pairs;
    // Cell -> indices into m_pairs, one per stacked tube segment. Empty = solid fill.
    std::unordered_map<CellId, std::vector<int>, CellIdHash> m_cell_pair_index;

    // Indexed by Layer::id().
    std::vector<LayerData> m_layer_data;

    // Cached zero-offset lattice (cell identity + topology); see topology_lattice().
    mutable std::unique_ptr<MagmaLattice> m_topology_lattice;

    InfillPattern         m_pattern = ipMagmaTriangle;
    const MagmaGeometry  *m_geometry = nullptr;

    // Config
    SpiralParams m_spiral_params;
    WindowSpec   m_window_spec;
    double m_cell_spacing;
    float  m_interior_width;
    float  m_line_width;
    double m_min_cap_clearance = 0.0; // min injection point -> opening boundary clearance
                                      // (injection nozzle flat radius)
    float  m_layer_height;            // nominal config layer height (fallback)
    float  m_min_layer_height;        // smallest layer height actually printed in this object
    double m_max_tube_height_mm;      // max tube height in mm (drives boundary placement)
    double m_min_tube_height_mm;      // min tube height in mm
    int    m_num_layers;
    // First m_layer_data index that is an object layer; non-zero with a raft, since
    // Layer::id() includes the raft layers.
    int    m_first_layer_id = 0;
    MagmaTubeSolverMode m_solver_mode = MagmaTubeSolverMode::CoverageStagger;
    double m_solver_timeout = 20.0;
    MagmaInjectionEdgePref m_injection_edge_pref = MagmaInjectionEdgePref::Interior;

    std::vector<std::string> m_warning_messages;  // slicing warnings to surface, in order found
    std::vector<int> m_injection_layer_ids;  // sorted, deduplicated cap layer IDs
    int    m_injection_extruder = 0;
};

using MagmaTubeMapPtr = std::unique_ptr<MagmaTubeMap>;

} // namespace magma
} // namespace Slic3r

#endif // slic3r_Magma_MagmaTubeMap_hpp_
