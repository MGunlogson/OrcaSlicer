#ifndef slic3r_Magma_MagmaTriangleCell_hpp_
#define slic3r_Magma_MagmaTriangleCell_hpp_

#include "../libslic3r.h"
#include "../Point.hpp"
#include "../BoundingBox.hpp"

#include <vector>
#include <array>
#include <cmath>
#include <functional>

namespace Slic3r {
namespace magma {

// ============================================================================
// Triangle Grid Constants
// ============================================================================

constexpr double SQRT3 = 1.7320508075688772935;
constexpr double INV_SQRT3 = 0.57735026918962576451;

// ============================================================================
// Triangle Grid Geometry Helpers
// ============================================================================

// Triangle side length (edge length) from cell spacing
// cell_spacing = distance between parallel lines
// side_length = cell_spacing * 2 / sqrt(3)
inline double triangle_side_length(double cell_spacing) {
    return cell_spacing * 2.0 * INV_SQRT3;
}

// Cell spacing from interior width and line width
// cell_spacing = interior_width + line_width (center-to-center distance)
inline double cell_spacing_from_geometry(double interior_width, double line_width) {
    return interior_width + line_width;
}

// Auto-calculate interior width from nozzle diameter
// Cell hole = nozzle + 0.2mm clearance for sealing during injection
double calculate_auto_interior_width(double nozzle_diameter);

// Auto-calculate window height from geometry
// Returns layers (minimum 2)
// Based on matching tube cross-section area to window area
int calculate_auto_window_height(double interior_width, double line_width, double layer_height);

// ============================================================================
// Triangle Cell Coordinate System
// ============================================================================

// Triangle cell coordinate system using (a, b, c) integer coordinates
//
// The triangular grid uses three axes at 60° angles, creating a dual coordinate system:
// - Up triangles: a + b + c == 2 (pointing up: △)
// - Down triangles: a + b + c == 1 (pointing down: ▽)
//
// This coordinate system allows:
// - Efficient neighbor lookups
// - Deterministic stagger level calculation: (a + b + c) % num_levels
// - Consistent cell identification across layers
//
// Geometry:
// - cell_spacing: distance between parallel lines (center-to-center)
// - Triangle side length = cell_spacing * 2 / sqrt(3)
// - Triangle height (altitude) = cell_spacing
//
struct TriangleCell {
    int a, b, c;  // Triangle coordinates

    TriangleCell() : a(0), b(0), c(0) {}
    TriangleCell(int a_, int b_, int c_) : a(a_), b(b_), c(c_) {}

    // Check if this is an upward-pointing triangle (△)
    bool is_up() const { return (a + b + c) == 2; }

    // Check if this is a downward-pointing triangle (▽)
    bool is_down() const { return (a + b + c) == 1; }

    // Get stagger level for window timing (0 to num_levels-1)
    int stagger_level(int num_levels) const {
        // Use all coordinates for better distribution
        int sum = a + b + c;
        // For up triangles sum=2, for down sum=1
        // Use a+b to differentiate within each type
        return ((a + b) % num_levels + num_levels) % num_levels;
    }

    // Get the paired cell for U-tube formation
    // Up triangle (a,b,c) pairs with down triangle (a+1,b,c) or (a,b+1,c) or (a,b,c+1)
    // Down triangle (a,b,c) pairs with up triangle (a-1,b,c) or (a,b-1,c) or (a,b,c-1)
    // The shared edge determines which pairing - default to the 'a' axis pairing
    TriangleCell get_paired_cell() const {
        if (is_up()) {
            // Up triangle (sum=2) pairs with down triangle (sum=1) at (a-1, b, c)
            return TriangleCell(a - 1, b, c);
        } else {
            // Down triangle (sum=1) pairs with up triangle (sum=2) at (a+1, b, c)
            return TriangleCell(a + 1, b, c);
        }
    }

    bool operator==(const TriangleCell& other) const {
        return a == other.a && b == other.b && c == other.c;
    }

    bool operator!=(const TriangleCell& other) const {
        return !(*this == other);
    }

    bool operator<(const TriangleCell& o) const {
        if (a != o.a) return a < o.a;
        if (b != o.b) return b < o.b;
        return c < o.c;
    }

    // Adjacent cells sharing an edge.
    // Up triangle neighbors: decrement one coordinate by 1 → down triangles (sum=1).
    // Down triangle neighbors: increment one coordinate by 1 → up triangles (sum=2).
    std::array<TriangleCell, 3> neighbors() const {
        if (is_up())
            // Up (sum=2) neighbors are down triangles (sum=1)
            return {{ {a-1,b,c}, {a,b-1,c}, {a,b,c-1} }};
        else
            // Down (sum=1) neighbors are up triangles (sum=2)
            return {{ {a+1,b,c}, {a,b+1,c}, {a,b,c+1} }};
    }
};

// Hash functor for TriangleCell, suitable for use in unordered containers.
struct TriangleCellHash {
    size_t operator()(const TriangleCell &c) const {
        // Combine all three coordinates with bit mixing
        size_t h = std::hash<int>()(c.a);
        h ^= std::hash<int>()(c.b) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= std::hash<int>()(c.c) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

// ============================================================================
// TriangleLattice - Unified Coordinate System for Triangle Grid
// ============================================================================
//
// The triangle grid uses a skewed coordinate system where:
// - Rows are spaced cell_spacing apart vertically
// - Columns are spaced edge_length apart, but shift by edge_length/2 per row
//
// Lattice coordinates (lx, ly) map to world coordinates (px, py) via:
//   py = ly * cell_spacing
//   px = lx * edge_length + ly * edge_length * 0.5
//
// This class bundles cell_spacing and spiral offset together, providing
// all coordinate transformations needed for cell detection, vertex finding,
// and window placement.

class TriangleLattice {
public:
    // Construct a lattice with given cell spacing and optional spiral offset
    explicit TriangleLattice(double cell_spacing, double offset_x = 0.0, double offset_y = 0.0)
        : m_cell_spacing(cell_spacing)
        , m_edge_length(triangle_side_length(cell_spacing))
        , m_offset_x(offset_x)
        , m_offset_y(offset_y)
    {}

    // Accessors
    double cell_spacing() const { return m_cell_spacing; }
    double edge_length() const { return m_edge_length; }
    double offset_x() const { return m_offset_x; }
    double offset_y() const { return m_offset_y; }

    // ========================================================================
    // Coordinate Transformations
    // ========================================================================

    // Convert lattice coordinates to world coordinates
    Vec2d to_world(double lx, double ly) const {
        return Vec2d(
            lx * m_edge_length + ly * m_edge_length * 0.5 + m_offset_x,
            ly * m_cell_spacing + m_offset_y
        );
    }

    // Convert world coordinates to lattice coordinates
    // Returns (lattice_x, lattice_y) in the unshifted reference frame
    std::pair<double, double> to_lattice(double px, double py) const {
        double adjusted_x = px - m_offset_x;
        double adjusted_y = py - m_offset_y;
        double ly = adjusted_y / m_cell_spacing;
        double lx = (adjusted_x - ly * m_edge_length * 0.5) / m_edge_length;
        return {lx, ly};
    }

    // ========================================================================
    // Cell Operations
    // ========================================================================

    // Get the triangle cell containing a world point
    TriangleCell cell_at(double px, double py) const {
        auto [lx, ly] = to_lattice(px, py);

        int col = static_cast<int>(std::floor(lx));
        int row = static_cast<int>(std::floor(ly));

        double fx = lx - col;
        double fy = ly - row;

        // fx + fy < 1 means UP triangle, >= 1 means DOWN triangle
        bool is_up = (fx + fy) < 1.0;

        int c = is_up ? (2 - col - row) : (1 - col - row);
        return TriangleCell(col, row, c);
    }

    // ========================================================================
    // Vertex Operations
    // ========================================================================

    // Get the nearest lattice vertex to a world point
    // Uses 4-candidate search due to skewed coordinate system
    Vec2d nearest_vertex(double px, double py) const {
        auto [lx, ly] = to_lattice(px, py);

        int lx0 = static_cast<int>(std::floor(lx));
        int lx1 = lx0 + 1;
        int ly0 = static_cast<int>(std::floor(ly));
        int ly1 = ly0 + 1;

        Vec2d best = to_world(lx0, ly0);
        double best_dist_sq = (px - best.x()) * (px - best.x()) + (py - best.y()) * (py - best.y());

        auto check = [&](int cx, int cy) {
            Vec2d v = to_world(cx, cy);
            double d = (px - v.x()) * (px - v.x()) + (py - v.y()) * (py - v.y());
            if (d < best_dist_sq) { best_dist_sq = d; best = v; }
        };
        check(lx1, ly0);
        check(lx0, ly1);
        check(lx1, ly1);

        return best;
    }

    // ========================================================================
    // Cell Geometry Operations
    // ========================================================================

    // Get the 3 corners of a cell in world coordinates
    std::array<Vec2d, 3> cell_corners(const TriangleCell& cell) const;

    // Get the center of a cell in world coordinates
    Vec2d cell_center(const TriangleCell& cell) const;

    // ========================================================================
    // Edge-Aware Generation Support
    // ========================================================================

    // Get the UP triangle cell for a horizontal edge at given row
    // Horizontal edges are bases of UP triangles, so this returns the cell that
    // "owns" this shared edge for window determination purposes.
    //
    // row: row index from horizontal line generation
    // x_center: X coordinate of line midpoint (to resolve which cell column)
    TriangleCell cell_for_horizontal_row(int row, double x_center) const;

    // Enumerate all unique cells within a bounding box.
    // Iterates lattice (col, row) coordinates covering the bbox,
    // generates both up and down triangles for each position.
    std::vector<TriangleCell> enumerate_cells(const BoundingBox& bbox) const;

private:
    double m_cell_spacing;
    double m_edge_length;
    double m_offset_x;
    double m_offset_y;
};

// Window specification for U-tube formation
// Windows are gaps in the shared walls between paired cells that create
// connected U-tubes. The stagger pattern prevents weak planes.
//
// Key relationships:
// - tube_height_layers = total height of one U-tube segment
// - window_height_layers = gap height at bottom of each tube segment
// - The "period" between windows equals tube_height_layers
// - Stagger offsets different cells' windows within the tube height cycle
//
// Example with tube_height=10, window_height=2, stagger_levels=3:
//   Stagger 0: windows at layers 0-1, 10-11, 20-21...
//   Stagger 1: windows at layers 3-4, 13-14, 23-24...
//   Stagger 2: windows at layers 6-7, 16-17, 26-27...
struct WindowSpec {
    int window_height_layers = 2;      // Gap height in layers
    int tube_height_layers = 10;       // Total tube segment height (= window period)
    int stagger_levels = 3;            // Number of Z-offset groups

    // Default constructor with defaults
    WindowSpec() = default;

    // Construct with explicit values
    WindowSpec(int window_height, int tube_height, int stagger)
        : window_height_layers(window_height)
        , tube_height_layers(tube_height)
        , stagger_levels(stagger)
    {}

    // Construct from config values with auto-calculation support
    // config_window_height = 0 means auto-calculate
    // config_tube_height_layers and config_tube_height_mm work together (max wins)
    static WindowSpec from_config(
        int config_window_height,           // 0 = auto-calculate
        int config_tube_height_layers,      // from magma_tube_height_layers
        float config_tube_height_mm,        // from magma_tube_height (0 = ignore)
        float interior_width,               // for auto window height calc
        float line_width,                   // for auto window height calc
        float layer_height,                 // for conversions
        int config_stagger_levels = 3       // number of stagger groups
    );

};

} // namespace magma
} // namespace Slic3r

#endif // slic3r_Magma_MagmaTriangleCell_hpp_
