#ifndef slic3r_Magma_MagmaHexCell_hpp_
#define slic3r_Magma_MagmaHexCell_hpp_

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

// Pure-hexagon (honeycomb) lattice of pointy-top hexagons: vertical left/right edges, four
// slanted edges, six neighbours. Adjacent hexes pair into a 2-cell U-tube with a window in their
// shared edge, as for triangle/rectilinear (no vents, unlike tri-hex).
//
// CellId carries (q, r, 0, 0): axial hex coordinates.
//
// The toolpath (FillMagmaHoneycomb) traces each vertical edge twice, once per adjacent column,
// so vertical walls are two beads thick. HexLattice folds that into its pitch so the open
// hexagon prints regular at flat-to-flat = interior.

// Hex edge length from flat-to-flat spacing: s = e*sqrt3.
inline double hex_edge_length(double cell_spacing) { return cell_spacing * INV_SQRT3; }

struct HexagonGeometry final : public MagmaGeometry
{
    // Insetting a 120deg corner by lw/2 shortens each edge by lw/sqrt3, so the open edge is
    // (spacing - lw)/sqrt3. Same expression as HexLattice::m_edge; keep them in agreement.
    double open_edge_length(double spacing, double line_width) const override {
        return hex_edge_length(spacing - line_width);
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

    // Honeycomb vertices are degree-3 line ends, not crossings, so there is no overlap.
    // The doubled vertical walls are captured by the measured cavity footprint.
    double vertex_overlap_excess_area(double /*line_width*/) const override {
        return 0.0;
    }

    // Open hex area / open shared-edge length.
    double auto_window_height(double interior_width, double line_width) const override {
        const double spacing   = interior_width + line_width;
        const double open_edge = open_edge_length(spacing, line_width);
        return open_edge > 0.0 ? inset_open_area(spacing, line_width) / open_edge : 0.0;
    }

    int max_neighbors() const override { return 6; }
};

inline const MagmaGeometry& hexagon_geometry() {
    static const HexagonGeometry s_geom;
    return s_geom;
}

// Pointy-top axial tiling: centre x = m_sx*(q + r/2) + ox, y = m_row*r + oy.
class HexLattice : public MagmaLattice {
public:
    HexLattice() = default;
    // The open hexagon has flat-to-flat = interior (cell_spacing - line_width); walls sit
    // outside it. Vertical walls are two beads (X pitch interior + 2*lw); slanted walls are
    // one bead, whose lw/2 offset lifts each apex by lw/sqrt3.
    explicit HexLattice(double cell_spacing, double offset_x = 0.0, double offset_y = 0.0,
                        double line_width = 0.0)
        : m_cell_spacing(cell_spacing)
        , m_edge(hex_edge_length(std::max(0.0, cell_spacing - std::max(0.0, line_width))))  // e = interior/√3
        , m_sx(cell_spacing + std::max(0.0, line_width))                                    // interior + 2·lw
        , m_vtop(m_edge + std::max(0.0, line_width) * INV_SQRT3)
        , m_row(1.5 * m_edge + std::max(0.0, line_width) * INV_SQRT3)
        , m_offset_x(offset_x)
        , m_offset_y(offset_y)
    {}

    // ---- topology: 6 edge-sharing neighbours (pointy-top axial) ----
    std::vector<CellId> neighbors(const CellId &c) const override {
        const int q = c.a, r = c.b;
        return { CellId(q+1, r,   0), CellId(q+1, r-1, 0), CellId(q,   r-1, 0),
                 CellId(q-1, r,   0), CellId(q-1, r+1, 0), CellId(q,   r+1, 0) };
    }
    bool is_up(const CellId &) const override { return false; }  // hexes have no parity
    int  max_neighbors() const override { return 6; }

    // ---- cell geometry ----
    Vec2d cell_center(const CellId &c) const override { return center(c.a, c.b); }

    std::vector<Vec2d> cell_corners(const CellId &c) const override {
        const Vec2d ctr = center(c.a, c.b);
        const double hx = m_sx * 0.5;              // half X pitch = interior/2 + lw
        const double hy = m_edge * 0.5;            // half open edge length
        // CCW from the upper-right corner (30deg).
        return { Vec2d(ctr.x() + hx, ctr.y() + hy),   //  30 deg (upper right)
                 Vec2d(ctr.x(),      ctr.y() + m_vtop),//  90 deg (top vertex, extended)
                 Vec2d(ctr.x() - hx, ctr.y() + hy),   // 150 deg (upper left)
                 Vec2d(ctr.x() - hx, ctr.y() - hy),   // 210 deg (lower left)
                 Vec2d(ctr.x(),      ctr.y() - m_vtop),// 270 deg (bottom vertex, extended)
                 Vec2d(ctr.x() + hx, ctr.y() - hy) }; // 330 deg (lower right)
    }

    CellId cell_at(double px, double py) const override {
        auto [fq, fr] = to_lattice(px, py);
        return axial_round(fq, fr);
    }

    std::vector<CellId> enumerate_cells(const BoundingBox &bbox) const override {
        double min_x = unscale<double>(bbox.min.x()), min_y = unscale<double>(bbox.min.y());
        double max_x = unscale<double>(bbox.max.x()), max_y = unscale<double>(bbox.max.y());
        double fq[4], fr[4];
        const double xs[4] = { min_x, max_x, min_x, max_x };
        const double ys[4] = { min_y, min_y, max_y, max_y };
        for (int k = 0; k < 4; ++k) { auto p = to_lattice(xs[k], ys[k]); fq[k] = p.first; fr[k] = p.second; }
        int q0 = int(std::floor(*std::min_element(fq, fq + 4))) - 1;
        int q1 = int(std::ceil (*std::max_element(fq, fq + 4))) + 1;
        int r0 = int(std::floor(*std::min_element(fr, fr + 4))) - 1;
        int r1 = int(std::ceil (*std::max_element(fr, fr + 4))) + 1;
        std::vector<CellId> cells;
        cells.reserve(size_t(std::max(0, q1 - q0 + 1)) * size_t(std::max(0, r1 - r0 + 1)));
        for (int r = r0; r <= r1; ++r)
            for (int q = q0; q <= q1; ++q)
                cells.emplace_back(q, r, 0);
        return cells;
    }

    // ---- coordinate transforms (pointy-top axial) ----
    Vec2d to_world(double lq, double lr) const override {
        return Vec2d(m_sx * (lq + lr * 0.5) + m_offset_x,
                     m_row * lr + m_offset_y);
    }
    std::pair<double, double> to_lattice(double px, double py) const override {
        double r = (py - m_offset_y) / m_row;
        double q = ((px - m_offset_x) / m_sx) - r * 0.5;
        return { q, r };
    }

    // ---- accessors ----
    double cell_spacing() const override { return m_cell_spacing; }
    double edge_length()  const override { return m_edge; }
    double offset_x()     const override { return m_offset_x; }
    double offset_y()     const override { return m_offset_y; }

private:
    Vec2d center(int q, int r) const { return to_world(double(q), double(r)); }

    // Round fractional axial (q, r) to the nearest hex via cube rounding.
    static CellId axial_round(double q, double r) {
        double x = q, z = r, y = -x - z;
        double rx = std::round(x), ry = std::round(y), rz = std::round(z);
        double dx = std::abs(rx - x), dy = std::abs(ry - y), dz = std::abs(rz - z);
        if (dx > dy && dx > dz)      rx = -ry - rz;
        else if (dy > dz)            ry = -rx - rz;
        else                         rz = -rx - ry;
        return CellId(int(rx), int(rz), 0);
    }

    double m_cell_spacing = 0.0;
    double m_edge = 0.0;        // open hex edge e = interior/sqrt3
    double m_sx   = 0.0;        // X pitch = interior + 2*line_width
    double m_vtop = 0.0;        // top/bottom vertex height = e + line_width/sqrt3
    double m_row  = 0.0;        // row spacing = 1.5*e + line_width/sqrt3
    double m_offset_x = 0.0;
    double m_offset_y = 0.0;
};

} // namespace magma
} // namespace Slic3r

#endif // slic3r_Magma_MagmaHexCell_hpp_
