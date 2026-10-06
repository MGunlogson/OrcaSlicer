#ifndef slic3r_Magma_MagmaRectilinearCell_hpp_
#define slic3r_Magma_MagmaRectilinearCell_hpp_

#include "../libslic3r.h"
#include "../Point.hpp"
#include "../BoundingBox.hpp"
#include "MagmaGeometry.hpp"
#include "MagmaLattice.hpp"
#include "MagmaCell.hpp"

#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>

namespace Slic3r {
namespace magma {

// Square lattice: two perpendicular line families cell_spacing apart. Walls are single beads on
// the grid lines, so the open square has side spacing - line_width. Seals better than the
// triangle (circumscribed/inscribed ratio sqrt2 vs 2) and prints two line families instead of 3.

constexpr double SQRT2 = 1.4142135623730951;

struct SquareGeometry final : public MagmaGeometry
{
    double open_edge_length(double spacing, double line_width) const override { return spacing - line_width; }

    // Open square of side spacing - line_width.
    double inset_open_area(double spacing, double line_width) const override {
        double s = spacing - line_width;
        return s > 0.0 ? s * s : 0.0;
    }

    // Circumscribed diameter of the open square: side * sqrt2.
    double opening_diameter(double spacing, double line_width) const override {
        double s = spacing - line_width;
        return s > 0.0 ? s * SQRT2 : 0.0;
    }

    double inscribed_radius(double interior_width, double) const override { return interior_width * 0.5; }

    // One 90deg crossing per cell, each double-depositing an lw x lw square.
    double vertex_overlap_excess_area(double line_width) const override {
        return line_width * line_width;
    }

    // Open area / open edge, which reduces to interior_width.
    double auto_window_height(double interior_width, double line_width) const override {
        const double spacing   = interior_width + line_width;
        const double open_edge = open_edge_length(spacing, line_width);
        return open_edge > 0.0 ? inset_open_area(spacing, line_width) / open_edge : 0.0;
    }

    int max_neighbors() const override { return 4; }
};

inline const MagmaGeometry& square_geometry() {
    static const SquareGeometry s_geom;
    return s_geom;
}

// CellId (a, b) = (column, row). Cell (a, b) is the unit square [a, a+1] x [b, b+1] in
// lattice space; world = lattice * cell_spacing + offset.
class RectilinearLattice : public MagmaLattice {
public:
    RectilinearLattice() = default;
    explicit RectilinearLattice(double cell_spacing, double offset_x = 0.0, double offset_y = 0.0)
        : m_cell_spacing(cell_spacing), m_offset_x(offset_x), m_offset_y(offset_y) {}

    // ---- topology ----
    std::vector<CellId> neighbors(const CellId &c) const override {
        return { CellId(c.a - 1, c.b, 0), CellId(c.a + 1, c.b, 0),
                 CellId(c.a, c.b - 1, 0), CellId(c.a, c.b + 1, 0) };
    }
    bool is_up(const CellId &) const override { return false; }  // squares have no parity
    int  max_neighbors() const override { return 4; }

    // ---- cell geometry ----
    CellId cell_at(double px, double py) const override {
        auto [lx, ly] = to_lattice(px, py);
        return CellId(int(std::floor(lx)), int(std::floor(ly)), 0);
    }
    std::vector<Vec2d> cell_corners(const CellId &c) const override {
        return { to_world(c.a,     c.b),
                 to_world(c.a + 1, c.b),
                 to_world(c.a + 1, c.b + 1),
                 to_world(c.a,     c.b + 1) };
    }
    Vec2d cell_center(const CellId &c) const override {
        return to_world(c.a + 0.5, c.b + 0.5);
    }
    std::vector<CellId> enumerate_cells(const BoundingBox &bbox) const override {
        std::vector<CellId> cells;
        double min_x = unscale<double>(bbox.min.x());
        double min_y = unscale<double>(bbox.min.y());
        double max_x = unscale<double>(bbox.max.x());
        double max_y = unscale<double>(bbox.max.y());
        auto [lx0, ly0] = to_lattice(min_x, min_y);
        auto [lx1, ly1] = to_lattice(max_x, max_y);
        int i0 = int(std::floor(std::min(lx0, lx1))) - 1;
        int i1 = int(std::ceil (std::max(lx0, lx1))) + 1;
        int j0 = int(std::floor(std::min(ly0, ly1))) - 1;
        int j1 = int(std::ceil (std::max(ly0, ly1))) + 1;
        cells.reserve(size_t(std::max(0, i1 - i0 + 1)) * size_t(std::max(0, j1 - j0 + 1)));
        for (int j = j0; j <= j1; ++j)
            for (int i = i0; i <= i1; ++i)
                cells.emplace_back(i, j, 0);
        return cells;
    }

    // ---- coordinate transforms ----
    Vec2d to_world(double lx, double ly) const override {
        return Vec2d(lx * m_cell_spacing + m_offset_x,
                     ly * m_cell_spacing + m_offset_y);
    }
    std::pair<double, double> to_lattice(double px, double py) const override {
        return { (px - m_offset_x) / m_cell_spacing,
                 (py - m_offset_y) / m_cell_spacing };
    }

    // ---- accessors ----
    double cell_spacing() const override { return m_cell_spacing; }
    double edge_length()  const override { return m_cell_spacing; }
    double offset_x()     const override { return m_offset_x; }
    double offset_y()     const override { return m_offset_y; }

private:
    double m_cell_spacing = 0.0;
    double m_offset_x = 0.0;
    double m_offset_y = 0.0;
};

} // namespace magma
} // namespace Slic3r

#endif // slic3r_Magma_MagmaRectilinearCell_hpp_
