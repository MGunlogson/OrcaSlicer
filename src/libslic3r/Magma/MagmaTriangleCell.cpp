#include "MagmaTriangleCell.hpp"

#include <cmath>
#include <algorithm>

namespace Slic3r {
namespace magma {

WindowSpec WindowSpec::from_config(
    const MagmaGeometry& geom,
    float config_window_height_mm,
    float interior_width,
    float line_width,
    float layer_height)
{
    WindowSpec spec;

    if (config_window_height_mm <= 0) {
        // Auto: the area-matched height scaled by MAGMA_AUTO_WINDOW_FRACTION, capped at the
        // window's width -- a window taller than it is wide lets plastic short-circuit across
        // its top instead of running down one leg and up the other. One layer is added so the
        // window spans at least a full printed layer.
        const double spacing      = double(interior_width) + double(line_width);
        const double window_width = geom.open_edge_length(spacing, double(line_width));
        double       h            = geom.auto_window_height(interior_width, line_width)
                                    * MAGMA_AUTO_WINDOW_FRACTION;
        if (window_width > 0.0)
            h = std::min(h, window_width);
        spec.window_height_mm = float(h) + std::max(0.0f, layer_height);
    } else {
        spec.window_height_mm = config_window_height_mm;
    }

    return spec;
}


// ============================================================================
// TriangleLattice Method Implementations
// ============================================================================

Vec2d TriangleLattice::cell_center(const CellId& cell) const
{
    // Centroid of the triangle's vertices (see cell_corners) in lattice space.
    if (is_up(cell)) {
        return to_world(cell.a + 1.0 / 3.0, cell.b + 1.0 / 3.0);
    } else {
        return to_world(cell.a + 2.0 / 3.0, cell.b + 2.0 / 3.0);
    }
}

std::vector<Vec2d> TriangleLattice::cell_corners(const CellId& cell) const
{
    // Up triangle (a,b): lattice (a,b), (a+1,b), (a,b+1).
    // Down triangle (a,b): lattice (a+1,b), (a,b+1), (a+1,b+1).
    if (is_up(cell)) {
        return { to_world(cell.a, cell.b),
                 to_world(cell.a + 1, cell.b),
                 to_world(cell.a, cell.b + 1) };
    } else {
        return { to_world(cell.a + 1, cell.b),
                 to_world(cell.a, cell.b + 1),
                 to_world(cell.a + 1, cell.b + 1) };
    }
}

std::vector<CellId> TriangleLattice::enumerate_cells(const BoundingBox& bbox) const
{
    std::vector<CellId> cells;

    double min_x = unscale<double>(bbox.min.x());
    double min_y = unscale<double>(bbox.min.y());
    double max_x = unscale<double>(bbox.max.x());
    double max_y = unscale<double>(bbox.max.y());

    int row_min = static_cast<int>(std::floor((min_y - m_offset_y) / m_cell_spacing)) - 1;
    int row_max = static_cast<int>(std::ceil((max_y - m_offset_y) / m_cell_spacing)) + 1;

    for (int row = row_min; row <= row_max; ++row) {
        // The lattice is skewed, so the column range depends on the row.
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
