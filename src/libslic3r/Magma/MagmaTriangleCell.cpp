#include "MagmaTriangleCell.hpp"

#include <cmath>
#include <algorithm>

namespace Slic3r {
namespace magma {

// Triangle grid geometry:
// - Lines are spaced cell_spacing apart
// - Triangle side length s = cell_spacing * 2 / sqrt(3)
// - Triangle height (altitude) h = cell_spacing
//
// Coordinate system basis vectors (for converting (a,b,c) to (x,y)):
// - a axis: (1, 0) scaled by s/2
// - b axis: (cos(60°), sin(60°)) = (0.5, sqrt(3)/2) scaled by s/2
// - c axis: (-cos(60°), sin(60°)) = (-0.5, sqrt(3)/2) scaled by s/2
//
// For a given cell_spacing:
// s = cell_spacing * 2 / sqrt(3)
// basis_a = (s/2, 0) = (cell_spacing / sqrt(3), 0)
// basis_b = (s/4, s*sqrt(3)/4) = (cell_spacing / (2*sqrt(3)), cell_spacing/2)
// basis_c = (-s/4, s*sqrt(3)/4) = (-cell_spacing / (2*sqrt(3)), cell_spacing/2)

// ============================================================================
// Geometry Helper Implementations
// ============================================================================

double calculate_auto_interior_width(double nozzle_diameter)
{
    return nozzle_diameter + 0.2;
}

int calculate_auto_window_height(double interior_width, double line_width, double min_layer_height)
{
    constexpr double triangle_factor = 0.433;  // sqrt(3)/4
    double window_height_mm = triangle_factor * interior_width * interior_width / line_width;
    // Use min_layer_height for conservative conversion: more layers = wider physical window.
    int layers = static_cast<int>(std::ceil(window_height_mm / min_layer_height)) + 1;
    return std::max(2, layers);  // Minimum 2 layers
}

// Calculate effective tube height from both config values (layers and mm).
// Uses min_layer_height for conservative mm→layers conversion.
static int calculate_tube_height_layers(int config_layers, float config_mm, float min_layer_height)
{
    int from_layers = config_layers;
    int from_mm = (config_mm > 0) ? static_cast<int>(std::ceil(config_mm / min_layer_height)) : 0;
    return std::max(from_layers, from_mm);
}

WindowSpec WindowSpec::from_config(
    int config_window_height,
    int config_tube_height_layers,
    float config_tube_height_mm,
    float interior_width,
    float line_width,
    float layer_height,
    int config_stagger_levels,
    float min_layer_height)
{
    WindowSpec spec;
    // Use min_layer_height if available and smaller, otherwise nominal
    float effective_lh = (min_layer_height > 0 && min_layer_height < layer_height)
        ? min_layer_height : layer_height;

    // Window height: auto-calculate if config is 0
    if (config_window_height <= 0) {
        spec.window_height_layers = calculate_auto_window_height(interior_width, line_width, effective_lh);
    } else {
        spec.window_height_layers = config_window_height;
    }

    // Tube height: auto-calculate if both configs are 0, otherwise take maximum of configured values
    spec.tube_height_layers = calculate_tube_height_layers(config_tube_height_layers, config_tube_height_mm, effective_lh);

    // Auto: derive from window height (enough solid wall above and below each window)
    int min_tube_height = spec.window_height_layers * 2 + 2;
    if (spec.tube_height_layers <= 0) {
        spec.tube_height_layers = min_tube_height;
    } else if (spec.tube_height_layers < min_tube_height) {
        // Configured value too small — clamp to structural minimum
        spec.tube_height_layers = min_tube_height;
    }

    // Stagger levels: limit to what fits in the tube height
    // Max stagger levels = tube_height / window_height (non-overlapping windows)
    int max_stagger = spec.tube_height_layers / spec.window_height_layers;
    spec.stagger_levels = std::clamp(config_stagger_levels, 1, std::max(1, max_stagger));

    return spec;
}


// ============================================================================
// TriangleLattice Method Implementations
// ============================================================================

Vec2d TriangleLattice::cell_center(const TriangleCell& cell) const
{
    // Centroid of the triangle within the lattice parallelogram cell (a, b).
    // Up triangle vertices: lattice (a,b), (a+1,b), (a,b+1) → centroid at (a+1/3, b+1/3)
    // Down triangle vertices: lattice (a+1,b), (a,b+1), (a+1,b+1) → centroid at (a+2/3, b+2/3)
    if (cell.is_up()) {
        return to_world(cell.a + 1.0 / 3.0, cell.b + 1.0 / 3.0);
    } else {
        return to_world(cell.a + 2.0 / 3.0, cell.b + 2.0 / 3.0);
    }
}

std::array<Vec2d, 3> TriangleLattice::cell_corners(const TriangleCell& cell) const
{
    // Use exact lattice vertex positions for correct geometry.
    // Up triangle at (a,b): vertices at lattice (a,b), (a+1,b), (a,b+1)
    // Down triangle at (a,b): vertices at lattice (a+1,b), (a,b+1), (a+1,b+1)
    if (cell.is_up()) {
        return {{ to_world(cell.a, cell.b),
                  to_world(cell.a + 1, cell.b),
                  to_world(cell.a, cell.b + 1) }};
    } else {
        return {{ to_world(cell.a + 1, cell.b),
                  to_world(cell.a, cell.b + 1),
                  to_world(cell.a + 1, cell.b + 1) }};
    }
}

TriangleCell TriangleLattice::cell_for_horizontal_row(int row, double x_center) const
{
    // Horizontal lines are the bases of UP triangles.
    // The Y position is: row * cell_spacing + offset_y
    // We nudge slightly upward to get the UP triangle above the line.

    double y_pos = row * m_cell_spacing + m_offset_y;
    double y_nudged = y_pos + m_cell_spacing * 0.01;  // Small nudge into the UP triangle

    TriangleCell cell = cell_at(x_center, y_nudged);

    // Ensure we return an UP triangle (shared edges belong to UP triangles by convention)
    if (!cell.is_up()) {
        cell = TriangleCell(cell.a, cell.b, 2 - cell.a - cell.b);
    }

    return cell;
}

std::vector<TriangleCell> TriangleLattice::enumerate_cells(const BoundingBox& bbox) const
{
    std::vector<TriangleCell> cells;

    double min_x = unscale<double>(bbox.min.x());
    double min_y = unscale<double>(bbox.min.y());
    double max_x = unscale<double>(bbox.max.x());
    double max_y = unscale<double>(bbox.max.y());

    // Row range is unskewed (ly = y / cell_spacing), so direct computation works
    int row_min = static_cast<int>(std::floor((min_y - m_offset_y) / m_cell_spacing)) - 1;
    int row_max = static_cast<int>(std::ceil((max_y - m_offset_y) / m_cell_spacing)) + 1;

    for (int row = row_min; row <= row_max; ++row) {
        // Column range depends on row due to lattice skew (lx = (x - ly*edge*0.5) / edge)
        double row_y = row * m_cell_spacing + m_offset_y;
        auto [lx_lo, _1] = to_lattice(min_x, row_y);
        auto [lx_hi, _2] = to_lattice(max_x, row_y);
        int col_min = static_cast<int>(std::floor(std::min(lx_lo, lx_hi))) - 1;
        int col_max = static_cast<int>(std::ceil(std::max(lx_lo, lx_hi))) + 1;

        for (int col = col_min; col <= col_max; ++col) {
            cells.emplace_back(col, row, 2 - col - row);  // up triangle
            cells.emplace_back(col, row, 1 - col - row);  // down triangle
        }
    }
    return cells;
}

} // namespace magma
} // namespace Slic3r
