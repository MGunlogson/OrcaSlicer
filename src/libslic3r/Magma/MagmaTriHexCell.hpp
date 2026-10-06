#ifndef slic3r_Magma_MagmaTriHexCell_hpp_
#define slic3r_Magma_MagmaTriHexCell_hpp_

#include "../libslic3r.h"
#include "../Point.hpp"
#include "../BoundingBox.hpp"
#include "MagmaGeometry.hpp"
#include "MagmaLattice.hpp"
#include "MagmaCell.hpp"
#include "MagmaTriangleCell.hpp"   // SQRT3, INV_SQRT3

#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <utility>

namespace Slic3r {
namespace magma {

// Trihexagonal (Kagome) lattice: the rectified triangular tiling. Hexagon cells sit at the
// vertices of an underlying triangular grid (edge g) and up/down triangle cells at the medial
// triangles of its faces. It is bipartite: a hex borders only triangles (6), a triangle only
// hexes (3). Hexes are the injection hubs, triangles the vents they feed.
// CellId carries (i, j, 0, kind) with kind = TriHexKind.
//
// The geometry below describes the hub. Per-cell open area comes from cell_corners() in the
// tube map, so there is no per-kind area method.

enum TriHexKind : uint8_t { THK_HEX = 0, THK_TRI_UP = 1, THK_TRI_DOWN = 2 };

// Floor for the auto window height (mm), which can come out sub-layer for small vents. A taller
// window cannot constrict flow; it only raises the minimum tube height. Not applied to a manual
// magma_window_height.
static constexpr double MAGMA_TRIHEX_MIN_WINDOW_MM = 1.0;

// Sizing (cell_spacing s = interior_width + line_width, interior = hex open flat-to-flat):
//   hex circumradius = trihex edge e  = s/sqrt3
//   underlying grid edge g            = 2e
//   vertex row y-step                 = g*sqrt3/2 = s
inline double trihex_edge_length(double cell_spacing) { return cell_spacing * INV_SQRT3; }   // e

struct HexGeometry final : public MagmaGeometry
{
    // Every window joins a hub hexagon to a vent triangle across one shared edge e. The vent's
    // 60deg corners shorten its side by lw*sqrt3 (the hub's 120deg corners only by lw/sqrt3),
    // so the vent side is the narrower end of the opening.
    double open_edge_length(double spacing, double line_width) const override {
        return trihex_edge_length(spacing) - line_width * SQRT3;
    }

    // Regular hexagon area 2*sqrt3*a^2, open apothem a = (s - lw)/2.
    double inset_open_area(double spacing, double line_width) const override {
        double a = (spacing - line_width) * 0.5;          // open apothem
        return a > 0.0 ? 2.0 * SQRT3 * a * a : 0.0;
    }

    // Circumscribed diameter of the open hexagon: 2 * (s - lw)/sqrt3.
    double opening_diameter(double spacing, double line_width) const override {
        double s = spacing - line_width;
        return s > 0.0 ? 2.0 * s * INV_SQRT3 : 0.0;
    }

    // Inradius of a regular hexagon is its apothem, interior/2.
    double inscribed_radius(double interior_width, double) const override { return interior_width * 0.5; }

    // Degree-4 vertices: two line families cross at 60deg, double-depositing a rhombus of
    // 2 w^2/sqrt3 per vertex. Charged per cell by corner count (hub 6 -> 1.5, vent 3 -> 0.75).
    double vertex_overlap_excess_area(double line_width) const override {
        return 2.0 * INV_SQRT3 * line_width * line_width;
    }

    // A window feeds only its vent (the hub is injected from the top), so the opening is
    // matched to the vent's open cross-section, not the much larger hub's. The vent is an
    // equilateral triangle of side e, inset by lw/2 to side e - lw*sqrt3.
    double auto_window_height(double interior_width, double line_width) const override {
        const double vent_side = open_edge_length(interior_width + line_width, line_width);
        if (vent_side <= 0.0)
            return 0.0;
        const double vent_area = (SQRT3 / 4.0) * vent_side * vent_side;
        return std::max(MAGMA_TRIHEX_MIN_WINDOW_MM, vent_area / vent_side);
    }

    int max_neighbors() const override { return 6; }   // hub's count; a vent has 3
};

inline const MagmaGeometry& trihex_geometry() {
    static const HexGeometry s_geom;
    return s_geom;
}

class TriHexLattice : public MagmaLattice {
public:
    TriHexLattice() = default;
    explicit TriHexLattice(double cell_spacing, double offset_x = 0.0, double offset_y = 0.0)
        : m_cell_spacing(cell_spacing)
        , m_edge(trihex_edge_length(cell_spacing))
        , m_offset_x(offset_x)
        , m_offset_y(offset_y)
    {}

    // ---- topology (bipartite hub<->vent) ----
    std::vector<CellId> neighbors(const CellId &c) const override {
        const int i = c.a, j = c.b;
        switch (c.kind) {
        case THK_HEX:
            return { {i,   j,   0, THK_TRI_UP},   {i-1, j,   0, THK_TRI_UP},   {i,   j-1, 0, THK_TRI_UP},
                     {i-1, j,   0, THK_TRI_DOWN}, {i,   j-1, 0, THK_TRI_DOWN}, {i-1, j-1, 0, THK_TRI_DOWN} };
        case THK_TRI_UP:    // up-tri(i,j) vertices V(i,j), V(i+1,j), V(i,j+1)
            return { {i,   j,   0, THK_HEX}, {i+1, j,   0, THK_HEX}, {i,   j+1, 0, THK_HEX} };
        case THK_TRI_DOWN:  // down-tri(i,j) vertices V(i+1,j), V(i,j+1), V(i+1,j+1)
        default:
            return { {i+1, j,   0, THK_HEX}, {i,   j+1, 0, THK_HEX}, {i+1, j+1, 0, THK_HEX} };
        }
    }

    bool is_up(const CellId &c) const override { return c.kind == THK_TRI_UP; }
    int  max_neighbors() const override { return 6; }

    // ---- cell geometry ----
    Vec2d cell_center(const CellId &c) const override {
        switch (c.kind) {
        case THK_HEX:
            return vertex(c.a, c.b);
        case THK_TRI_UP:
            return centroid3(vertex(c.a, c.b), vertex(c.a + 1, c.b), vertex(c.a, c.b + 1));
        case THK_TRI_DOWN:
        default:
            return centroid3(vertex(c.a + 1, c.b), vertex(c.a, c.b + 1), vertex(c.a + 1, c.b + 1));
        }
    }

    std::vector<Vec2d> cell_corners(const CellId &c) const override {
        if (c.kind == THK_HEX) {
            // 6 edge-midpoints to the grid-neighbours, CCW from 0deg.
            const int i = c.a, j = c.b;
            Vec2d v = vertex(i, j);
            const std::array<std::pair<int,int>,6> nb = {{
                {i+1, j}, {i, j+1}, {i-1, j+1}, {i-1, j}, {i, j-1}, {i+1, j-1} }};
            std::vector<Vec2d> out;
            out.reserve(6);
            for (auto &n : nb) out.push_back(0.5 * (v + vertex(n.first, n.second)));
            return out;
        }
        // Triangle vent = medial triangle (edge midpoints) of the underlying face.
        Vec2d A, B, C;
        if (c.kind == THK_TRI_UP) {
            A = vertex(c.a, c.b); B = vertex(c.a + 1, c.b); C = vertex(c.a, c.b + 1);
        } else {
            A = vertex(c.a + 1, c.b); B = vertex(c.a, c.b + 1); C = vertex(c.a + 1, c.b + 1);
        }
        return { 0.5 * (A + B), 0.5 * (B + C), 0.5 * (C + A) };
    }

    CellId cell_at(double px, double py) const override {
        // Find the underlying triangle and barycentric coords. A coord > 0.5 puts the point
        // in that vertex's hexagon; otherwise it is in the medial triangle.
        auto [lx, ly] = to_lattice(px, py);
        int col = int(std::floor(lx)), row = int(std::floor(ly));
        double fx = lx - col, fy = ly - row;
        if (fx + fy < 1.0) {
            // up-tri(col,row): bary (1-fx-fy, fx, fy) for V(col,row),V(col+1,row),V(col,row+1)
            double w0 = 1.0 - fx - fy, w1 = fx, w2 = fy;
            if (w0 > 0.5) return { col,   row,   0, THK_HEX };
            if (w1 > 0.5) return { col+1, row,   0, THK_HEX };
            if (w2 > 0.5) return { col,   row+1, 0, THK_HEX };
            return { col, row, 0, THK_TRI_UP };
        } else {
            // down-tri(col,row): bary for V(col+1,row),V(col,row+1),V(col+1,row+1)
            double w0 = 1.0 - fy, w1 = 1.0 - fx, w2 = fx + fy - 1.0;
            if (w0 > 0.5) return { col+1, row,   0, THK_HEX };
            if (w1 > 0.5) return { col,   row+1, 0, THK_HEX };
            if (w2 > 0.5) return { col+1, row+1, 0, THK_HEX };
            return { col, row, 0, THK_TRI_DOWN };
        }
    }

    std::vector<CellId> enumerate_cells(const BoundingBox &bbox) const override {
        double min_x = unscale<double>(bbox.min.x()), min_y = unscale<double>(bbox.min.y());
        double max_x = unscale<double>(bbox.max.x()), max_y = unscale<double>(bbox.max.y());
        // Map the four bbox corners to lattice space and take the spanning (i,j) box.
        double lxs[4], lys[4];
        const double xs[4] = { min_x, max_x, min_x, max_x };
        const double ys[4] = { min_y, min_y, max_y, max_y };
        for (int k = 0; k < 4; ++k) { auto p = to_lattice(xs[k], ys[k]); lxs[k] = p.first; lys[k] = p.second; }
        int i0 = int(std::floor(*std::min_element(lxs, lxs + 4))) - 1;
        int i1 = int(std::ceil (*std::max_element(lxs, lxs + 4))) + 1;
        int j0 = int(std::floor(*std::min_element(lys, lys + 4))) - 1;
        int j1 = int(std::ceil (*std::max_element(lys, lys + 4))) + 1;
        std::vector<CellId> cells;
        cells.reserve(size_t(std::max(0, i1 - i0 + 1)) * size_t(std::max(0, j1 - j0 + 1)) * 3);
        for (int j = j0; j <= j1; ++j)
            for (int i = i0; i <= i1; ++i) {
                cells.emplace_back(i, j, 0, THK_HEX);
                cells.emplace_back(i, j, 0, THK_TRI_UP);
                cells.emplace_back(i, j, 0, THK_TRI_DOWN);
            }
        return cells;
    }

    // ---- coordinate transforms ----
    // Underlying grid basis u1 = (g, 0), u2 = (g/2, s); vertex(i,j) = i*u1 + j*u2.
    Vec2d to_world(double lx, double ly) const override {
        double g = 2.0 * m_edge;
        return Vec2d(lx * g + ly * (g * 0.5) + m_offset_x,
                     ly * m_cell_spacing + m_offset_y);
    }
    std::pair<double, double> to_lattice(double px, double py) const override {
        double g = 2.0 * m_edge;
        double ly = (py - m_offset_y) / m_cell_spacing;
        double lx = ((px - m_offset_x) - ly * (g * 0.5)) / g;
        return { lx, ly };
    }

    // ---- accessors ----
    double cell_spacing() const override { return m_cell_spacing; }
    double edge_length()  const override { return m_edge; }
    double offset_x()     const override { return m_offset_x; }
    double offset_y()     const override { return m_offset_y; }

private:
    Vec2d vertex(int i, int j) const { return to_world(double(i), double(j)); }
    static Vec2d centroid3(const Vec2d &a, const Vec2d &b, const Vec2d &c) {
        return (a + b + c) / 3.0;
    }

    double m_cell_spacing = 0.0;
    double m_edge = 0.0;       // trihex wall length e = spacing/sqrt3
    double m_offset_x = 0.0;
    double m_offset_y = 0.0;
};

} // namespace magma
} // namespace Slic3r

#endif // slic3r_Magma_MagmaTriHexCell_hpp_
