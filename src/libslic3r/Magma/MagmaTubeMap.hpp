#ifndef slic3r_Magma_MagmaTubeMap_hpp_
#define slic3r_Magma_MagmaTubeMap_hpp_

#include "MagmaTriangleCell.hpp"
#include "MagmaSpiralOffset.hpp"
#include "../PrintConfig.hpp"

#include <map>
#include <memory>
#include <unordered_map>
#include <vector>
#include <climits>

namespace Slic3r {

class Layer;

namespace magma {

// Per-cell layer presence and interior area
struct CellPresence {
    int first_layer = INT_MAX;   // first layer where cell exists in magma region
    int last_layer  = -1;        // last layer where cell exists in magma region
    std::vector<bool>   layers;  // indexed [layer_id - first_layer]
    std::vector<double> areas;   // parallel, interior area in scaled^2 units

    bool   present(int layer_id) const;
    double area(int layer_id) const;
    void   mark_present(int layer_id, double area_val);
};

// A U-tube pair assignment
struct UTubePair {
    TriangleCell cell_a;          // injection side (up triangle)
    TriangleCell cell_b;          // vent side (down triangle)
    int  pair_start_layer;        // first layer of this tube segment
    int  pair_end_layer;          // last layer of this tube segment (inclusive)
    bool is_salvaged;             // true = greedy salvage, false = default pair
    double volume_mm3;            // injection volume (both cells combined)
};

// Which edge two adjacent triangle cells share.
enum class SharedEdge { Horizontal, Col60, Diag120 };

// Determine shared edge type between two adjacent cells.
// Cells differ in exactly one coordinate (a, b, or c).
inline SharedEdge shared_edge(const TriangleCell &a, const TriangleCell &b) {
    if (a.a != b.a) return SharedEdge::Col60;       // differ in a → 60° edge
    if (a.b != b.b) return SharedEdge::Horizontal;   // differ in b → horizontal
    return SharedEdge::Diag120;                       // differ in c → 120° edge
}

// All window gap data for a layer, returned by MagmaTubeMap::window_gaps().
//
// Gaps are organized by line family (horizontal, 60°, 120°). Each map key is
// the line index (row, column, or diagonal), and the value is a sorted,
// merged list of world-coordinate intervals where lines should be interrupted.
//
// - Horizontal: key = row b, intervals are X ranges
// - Col60:      key = column a, intervals are Y ranges
// - Diag120:    key = diagonal s (= a+b+1 of UP cell), intervals are Y ranges
struct WindowGaps {
    double cell_spacing;
    double offset_x, offset_y;  // spiral offset for this layer (mm)

    std::map<int, std::vector<std::pair<double, double>>> horiz;
    std::map<int, std::vector<std::pair<double, double>>> col60;
    std::map<int, std::vector<std::pair<double, double>>> diag120;

    bool empty() const { return horiz.empty() && col60.empty() && diag120.empty(); }
};

class MagmaTubeMap {
public:
    // Build the tube map. Called from PrintObject::infill() before TBB fill loop.
    static std::unique_ptr<MagmaTubeMap> build(
        const std::vector<Layer*> &layers,
        const PrintRegionConfig &config,
        const PrintConfig &print_config,
        const PrintObjectConfig &obj_config);

    // === Query interface (thread-safe, const) ===

    // Pre-computed window gap intervals for a layer. Returns world-space
    // intervals for all three line families (horizontal, 60°, 120°) where
    // infill lines should be interrupted. Shared-edge detection ensures
    // gaps appear on the correct edge between paired cells.
    WindowGaps window_gaps(int layer_id) const;

    // Is window open for this cell on this layer?
    bool is_window_open(const TriangleCell &cell, int layer_id) const;

    // Is this cell assigned to a tube (true) or should be solid fill (false)?
    bool is_paired(const TriangleCell &cell) const;

    // Accessors for shared params (so FillMagma uses identical values)
    const SpiralParams& spiral_params() const { return m_spiral_params; }
    double cell_spacing() const { return m_cell_spacing; }
    float  interior_width() const { return m_interior_width; }

    // For injection G-code (future)
    const std::vector<UTubePair>& u_tube_pairs() const { return m_pairs; }

    // Statistics for debug logging
    int num_cells() const { return static_cast<int>(m_cells.size()); }
    int num_default_pairs() const;
    int num_salvaged_pairs() const;
    int num_solid_cells() const;

private:
    MagmaTubeMap() = default;

    // Build phases
    void scan_layers(const std::vector<Layer*> &layers);
    void detect_constrictions();
    void assign_default_tubes();
    void assign_salvage_tubes();
    void compute_volumes(const std::vector<Layer*> &layers);

    // Data
    std::unordered_map<TriangleCell, CellPresence, TriangleCellHash> m_cells;
    std::vector<UTubePair> m_pairs;
    // Maps each cell to its tube pair indices. A cell with multiple stacked
    // segments (from phase-shifted boundaries) will have multiple entries.
    // Empty vector = solid fill (no tube partner found).
    // Absent = not yet processed by assign_default/salvage_tubes.
    std::unordered_map<TriangleCell, std::vector<int>, TriangleCellHash> m_cell_pair_index;

    // Config
    SpiralParams m_spiral_params;
    WindowSpec   m_window_spec;
    double m_cell_spacing;
    float  m_interior_width;
    float  m_line_width;
    float  m_layer_height;
    int    m_min_tube_height_layers;
    int    m_max_tube_height_layers;
    int    m_num_layers;
    bool   m_dual_infill_enabled;
};

using MagmaTubeMapPtr = std::unique_ptr<MagmaTubeMap>;

} // namespace magma
} // namespace Slic3r

#endif // slic3r_Magma_MagmaTubeMap_hpp_
