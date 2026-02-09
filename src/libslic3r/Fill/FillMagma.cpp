#include "FillMagma.hpp"
#include "FillBase.hpp"
#include "../ClipperUtils.hpp"
#include "../Polyline.hpp"
#include "../Magma/MagmaTubeMap.hpp"
#include "../Magma/MagmaSpiralOffset.hpp"

#include <cmath>
#include <algorithm>
#include <boost/log/trivial.hpp>

namespace Slic3r {

// ============================================================================
// Infill Direction — fixed angle, world-origin reference
// ============================================================================

std::pair<float, Point> FillMagmaTriangle::_infill_direction(const Surface *surface) const
{
    float out_angle = 0.f;
    if (surface->bridge_angle >= 0)
        out_angle = float(surface->bridge_angle);
    return std::make_pair(out_angle, Point(0, 0));
}

// ============================================================================
// Helpers
// ============================================================================

// Create a 2-point polyline from world coordinates (mm) for X, scaled Y.
static Polyline make_horiz_segment(double x0_mm, double x1_mm, coord_t y_scaled)
{
    Polyline pl;
    pl.points.push_back(Point(coord_t(scale_(x0_mm)), y_scaled));
    pl.points.push_back(Point(coord_t(scale_(x1_mm)), y_scaled));
    return pl;
}

// Subtract sorted, merged gap intervals from a range [lo, hi].
// Emits kept segments to `out` via the callback.
template<typename EmitFn>
static void subtract_gaps(double lo, double hi,
                          const std::vector<std::pair<double, double>> &gaps,
                          EmitFn emit)
{
    double cursor = lo;
    for (const auto &gap : gaps) {
        double gl = std::max(gap.first, lo);
        double gr = std::min(gap.second, hi);
        if (gr <= gl)
            continue;
        if (gl > cursor + 0.01)
            emit(cursor, gl);
        cursor = std::max(cursor, gr);
    }
    if (hi > cursor + 0.01)
        emit(cursor, hi);
}

// Split a line segment p0→p1 by removing Y intervals where gaps exist.
// Gaps are in world-space Y coordinates. The line must have non-zero Y span.
static void split_line_by_y_gaps(
    const Vec2d &p0, const Vec2d &p1,
    const std::vector<std::pair<double, double>> &y_gaps,
    Polylines &out)
{
    double dy = p1.y() - p0.y();
    if (std::abs(dy) < 1e-6) {
        // Degenerate — no Y span, can't split by Y
        Polyline pl;
        pl.points.push_back(Point(scale_(p0.x()), scale_(p0.y())));
        pl.points.push_back(Point(scale_(p1.x()), scale_(p1.y())));
        out.push_back(std::move(pl));
        return;
    }

    double y_lo = std::min(p0.y(), p1.y());
    double y_hi = std::max(p0.y(), p1.y());

    // Convert Y intervals to parametric t values along p0→p1
    // P(t) = p0 + t*(p1 - p0), t ∈ [0, 1]
    // t = (y - p0.y) / dy
    auto y_to_t = [&](double y) { return (y - p0.y()) / dy; };
    auto t_to_point = [&](double t) -> Vec2d {
        return Vec2d(p0.x() + t * (p1.x() - p0.x()),
                     p0.y() + t * dy);
    };

    // Collect gap t-intervals
    std::vector<std::pair<double, double>> t_gaps;
    for (const auto &gap : y_gaps) {
        double gl = std::max(gap.first, y_lo);
        double gr = std::min(gap.second, y_hi);
        if (gr <= gl)
            continue;
        double t0 = y_to_t(gl);
        double t1 = y_to_t(gr);
        if (t0 > t1) std::swap(t0, t1);
        t_gaps.push_back({t0, t1});
    }

    if (t_gaps.empty()) {
        Polyline pl;
        pl.points.push_back(Point(scale_(p0.x()), scale_(p0.y())));
        pl.points.push_back(Point(scale_(p1.x()), scale_(p1.y())));
        out.push_back(std::move(pl));
        return;
    }

    std::sort(t_gaps.begin(), t_gaps.end());

    subtract_gaps(0.0, 1.0, t_gaps, [&](double t_lo, double t_hi) {
        Vec2d a = t_to_point(t_lo);
        Vec2d b = t_to_point(t_hi);
        Polyline pl;
        pl.points.push_back(Point(scale_(a.x()), scale_(a.y())));
        pl.points.push_back(Point(scale_(b.x()), scale_(b.y())));
        out.push_back(std::move(pl));
    });
}

// ============================================================================
// Main Fill — Direct Lattice Generation
// ============================================================================

void FillMagmaTriangle::_fill_surface_single(
    const FillParams &params,
    unsigned int       thickness_layers,
    const std::pair<float, Point> &direction,
    ExPolygon          expolygon,
    Polylines         &polylines_out)
{
    if (!this->tube_map) {
        BOOST_LOG_TRIVIAL(error) << "FillMagmaTriangle: null tube_map on layer " << this->layer_id;
        return;
    }

    const double cs   = this->tube_map->cell_spacing();
    const int    layer = static_cast<int>(this->layer_id);

    // Build lattice with spiral offset for this layer
    magma::TriangleLattice lattice = magma::lattice_for_layer(
        cs, this->tube_map->spiral_params(), layer);

    const double edge  = lattice.edge_length();
    const double off_x = lattice.offset_x();
    const double off_y = lattice.offset_y();

    // Get window gaps (pre-computed by tube map, correct shared-edge detection)
    magma::WindowGaps gaps = this->tube_map->window_gaps(layer);

    // ---- Bounding box → lattice ranges ----

    BoundingBox bbox = expolygon.contour.bounding_box();
    double x_min = unscale<double>(bbox.min.x()) - edge;
    double x_max = unscale<double>(bbox.max.x()) + edge;
    double y_min = unscale<double>(bbox.min.y()) - cs;
    double y_max = unscale<double>(bbox.max.y()) + cs;

    // Row range (horizontal lines, constant b)
    int row_min = static_cast<int>(std::floor((y_min - off_y) / cs));
    int row_max = static_cast<int>(std::ceil((y_max - off_y) / cs));

    // Column range (60° lines, constant a)
    // to_world(a, b).x = a*edge + b*edge/2 + off_x
    // Conservative: use extreme b values that push x furthest
    int col_min = static_cast<int>(std::floor(
        (x_min - off_x - std::max(row_min, row_max) * edge * 0.5) / edge)) - 1;
    int col_max = static_cast<int>(std::ceil(
        (x_max - off_x - std::min(row_min, row_max) * edge * 0.5) / edge)) + 1;

    // Diagonal range (120° lines, constant s = a + b)
    int diag_min = col_min + row_min - 1;
    int diag_max = col_max + row_max + 1;

    // ---- Generate lines with gaps built in ----

    Polylines all_lines;

    // --- Horizontal lines (one per row b) ---
    for (int b = row_min; b <= row_max; ++b) {
        coord_t y_s = coord_t(scale_(b * cs + off_y));

        auto it = gaps.horiz.find(b);
        if (it == gaps.horiz.end()) {
            // No gaps — single spanning line
            all_lines.push_back(make_horiz_segment(x_min, x_max, y_s));
        } else {
            // Subtract gap X intervals from [x_min, x_max]
            subtract_gaps(x_min, x_max, it->second, [&](double lo, double hi) {
                all_lines.push_back(make_horiz_segment(lo, hi, y_s));
            });
        }
    }

    // --- 60° lines (one per column a) ---
    for (int a = col_min; a <= col_max; ++a) {
        Vec2d p0 = lattice.to_world(a, row_min - 1);
        Vec2d p1 = lattice.to_world(a, row_max + 1);

        auto it = gaps.col60.find(a);
        if (it == gaps.col60.end()) {
            Polyline pl;
            pl.points.push_back(Point(scale_(p0.x()), scale_(p0.y())));
            pl.points.push_back(Point(scale_(p1.x()), scale_(p1.y())));
            all_lines.push_back(std::move(pl));
        } else {
            split_line_by_y_gaps(p0, p1, it->second, all_lines);
        }
    }

    // --- 120° lines (one per diagonal s = a + b) ---
    for (int s = diag_min; s <= diag_max; ++s) {
        // Line goes through lattice vertices (s - b, b) as b increases.
        // Direction in world: to_world(s-b, b) → to_world(s-b-1, b+1)
        //   Δworld = (-edge + edge/2, cs) = (-edge/2, cs)
        Vec2d p0 = lattice.to_world(s - (row_min - 1), row_min - 1);
        Vec2d p1 = lattice.to_world(s - (row_max + 1), row_max + 1);

        auto it = gaps.diag120.find(s);
        if (it == gaps.diag120.end()) {
            Polyline pl;
            pl.points.push_back(Point(scale_(p0.x()), scale_(p0.y())));
            pl.points.push_back(Point(scale_(p1.x()), scale_(p1.y())));
            all_lines.push_back(std::move(pl));
        } else {
            split_line_by_y_gaps(p0, p1, it->second, all_lines);
        }
    }

    // ---- Clip to polygon and anchor ----

    // Clip all lines to the fill polygon boundary
    all_lines = intersection_pl(std::move(all_lines), expolygon);

    // Anchor to perimeter walls + optimize travel path
    chain_or_connect_infill(std::move(all_lines), expolygon,
                            polylines_out, this->spacing, params);
}

} // namespace Slic3r
