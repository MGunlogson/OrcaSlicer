#ifndef slic3r_Magma_MagmaTriangleCell_hpp_
#define slic3r_Magma_MagmaTriangleCell_hpp_

#include "../libslic3r.h"
#include "../Point.hpp"
#include "../BoundingBox.hpp"
#include "MagmaGeometry.hpp"
#include "MagmaCell.hpp"
#include "MagmaLattice.hpp"

#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <functional>
#include <utility>

namespace Slic3r {
namespace magma {

constexpr double SQRT3 = 1.7320508075688772935;
constexpr double INV_SQRT3 = 0.57735026918962576451;

// Side of the equilateral triangle whose altitude is cell_spacing (the line spacing).
inline double triangle_side_length(double cell_spacing) {
    return cell_spacing * 2.0 * INV_SQRT3;
}

// Center-to-center line spacing.
inline double cell_spacing_from_geometry(double interior_width, double line_width) {
    return interior_width + line_width;
}

// Open tube cross-section (mm^2): the triangle inset by half a bead on each edge.
inline double inset_triangle_area(double cell_spacing, double line_width) {
    double side = triangle_side_length(cell_spacing);
    double inset_side = side - line_width * SQRT3;
    return (inset_side > 0) ? (SQRT3 / 4.0) * inset_side * inset_side : 0.0;
}

// Injection seal geometry, shared by the injection g-code and Print::validate() so the predicted
// seal matches the emitted one. The nozzle seals a tube when its flat tip, plus the cone above it
// once pressed down, covers the opening.

constexpr double MAGMA_DEG2RAD           = 0.017453292519943295;
constexpr double MAGMA_SLAM_CLAMP        = 3.5;  // ceiling on seal + plunge depth (mm)

// Cone-diameter-to-cell-pitch ratio above which the nozzle crushes neighbouring cells. A geometry
// sanity bound, not a quality threshold: at 2.0 the cone spans a whole neighbouring cell. Observed
// wall damage tracks tube width, not this ratio.
constexpr double MAGMA_PITCH_ABSURD_RATIO = 2.0;

// Seal depth (mm) below which the seal is reported as unreliable: a shallow engagement has no
// tolerance left for layer-height and squish variation on the printed rim. Test prints leaked at
// 0.29 and held at 0.56.
constexpr double MAGMA_SEAL_DEPTH_MIN     = 0.40;

// Longest an injection may run (s). A sealed nozzle is a heat source against thin walls; held
// longer, it softens them, the lattice deforms and the seal leaks. Duration matters, not volume
// or tube size. The rate is capped by the filament's max volumetric speed, so the way to shorten
// an injection is less volume per tube; splitting it into bursts fails because the melt freezes
// into a plug.
constexpr double MAGMA_MAX_INJECTION_SECONDS = 1.5;

// Estimated seconds for one U-tube injection (a pair of cells). Ignores the window cavity and
// edge clipping; the slicer measures each tube's real cavity from the toolpath. Returns 0 when
// the rate is unknown so callers can skip the check.
inline double injection_seconds(double open_area_mm2, double tube_height_mm, double fill_factor,
                                double max_volumetric_mm3s) {
    if (max_volumetric_mm3s <= 0.0 || open_area_mm2 <= 0.0 || tube_height_mm <= 0.0)
        return 0.0;
    return (2.0 * open_area_mm2 * tube_height_mm * fill_factor) / max_volumetric_mm3s;
}
// Default magma_seal_press (mm): how far the nozzle descends past first contact with the cell.
// See corner_grip().
constexpr double MAGMA_SEAL_PRESS         = 0.1;
// Clearance added to the auto crater-iron start margin, beyond the r_flat needed to put the
// nozzle flat outside the crater mouth. Covers the displaced material that rode up the bevel
// and piled into a rim just outside the mouth.
constexpr double MAGMA_IRON_RIM_CLEARANCE = 0.3;
// Auto crater-iron radial step per spiral turn, as a fraction of the nozzle flat. The flat
// is also the width of the bevel doing the plowing, so the step must scale with it.
constexpr double MAGMA_IRON_STEP_FRACTION = 0.5;

// Flat-from-bore estimate, used only as the suggested starting value in the unmeasured-flat
// error; Print::validate() refuses to slice with the flat unset. Deliberately below the measured
// 2.83x of an E3D 0.6mm nozzle: overestimating sizes tubes the nozzle cannot seal, while
// underestimating only gives smaller tubes.
constexpr double MAGMA_FLAT_BORE_MULTIPLE = 2.5;
inline double resolve_nozzle_flat(double configured_flat, double nozzle_diameter) {
    return configured_flat > 0.0 ? configured_flat
                                 : MAGMA_FLAT_BORE_MULTIPLE * nozzle_diameter;
}

// Depth the cone must descend so it widens from `flat` to cover `opening_dia`.
inline double seal_depth_for_opening(double opening_dia, double flat, double cone_half_angle_deg) {
    double tan_t = std::tan(cone_half_angle_deg * MAGMA_DEG2RAD);
    double gap = opening_dia - flat;
    return (gap > 0.0 && tan_t > 1e-6) ? gap / (2.0 * tan_t) : 0.0;
}

// Diameter of the nozzle cone at a given depth below first contact with the print surface.
inline double cone_diameter_at(double depth, double flat, double cone_half_angle_deg) {
    return flat + 2.0 * std::max(0.0, depth) * std::tan(cone_half_angle_deg * MAGMA_DEG2RAD);
}

// Radial grip on the cell corners at the end of the injection, which is what holds the seal.
// The corners are the last part of the opening the cone covers, so at seal depth they are
// gripped only by the press; the plunge adds more as it descends. Independent of tube size.
inline double corner_grip(double seal_press, double plunge_depth, double cone_half_angle_deg) {
    return (std::max(0.0, seal_press) + std::max(0.0, plunge_depth))
         * std::tan(cone_half_angle_deg * MAGMA_DEG2RAD);
}

// Depth at which the cone covers the whole opening, plus `seal_press`. Reached by one fast Z move
// before any flow; the plunge goes deeper during the fill. `opening_dia` is the circumscribed
// circle, since the corners are the last thing the cone reaches. The depth follows from the tube
// width and the nozzle; it is reported, not configured.
inline double auto_seal_depth(double opening_dia, double flat, double cone_half_angle_deg,
                              double seal_press) {
    // First contact: zero when the flat already spans the opening.
    const double first_contact = seal_depth_for_opening(opening_dia, flat, cone_half_angle_deg);
    // At first contact the nozzle rests rather than grips; the press supplies the grip.
    return first_contact + std::max(0.0, seal_press);
}

// The plunge descends during the injection, on top of the seal depth; only MAGMA_SLAM_CLAMP
// bounds the sum. Print::validate() checks the result against the lattice pitch.
// Plunge the user asked for: the configured depth when plunging is on, else none.
inline double requested_plunge_depth(bool plunge_enabled, double plunge_depth_cfg) {
    return plunge_enabled ? std::max(0.0, plunge_depth_cfg) : 0.0;
}

inline double clamp_plunge_depth(double seal_depth, double plunge_depth) {
    return std::min(std::max(0.0, plunge_depth),
                    std::max(0.0, MAGMA_SLAM_CLAMP - std::max(0.0, seal_depth)));
}

// Crater-iron pass. MagmaInjection.cpp emits it and Print.cpp's heat-spread timing model clocks
// it, so both take the geometry from here.
struct CraterIron
{
    double start_R = 0.0;  // spiral start radius: crater footprint plus margin
    int    turns   = 1;    // spiral turns from start_R in to the centre

    // Path length: one rim revolution, an Archimedean spiral to the centre (mean radius
    // start_R/2 per turn), and the flatten stroke out, across and back (2 flat widths).
    double length(double flat) const { return 2.0 * PI * start_R + PI * turns * start_R + 2.0 * flat; }
};

// `press_depth` is how far below the surface the nozzle went (seal + plunge). The crater
// footprint is where the cone reached at that depth. Margin 0 = auto: the spiral plows the
// displaced rim with the bevel, so the flat has to start entirely outside it (margin >= flat
// radius) plus room for the piled material. Turns 0 = auto: the radial step scales with the
// flat, because the flat is the width of the bevel doing the plowing.
inline CraterIron plan_crater_iron(double flat, double cone_half_angle_deg, double press_depth,
                                   double margin_cfg, int turns_cfg)
{
    const double crater_r = 0.5 * cone_diameter_at(press_depth, flat, cone_half_angle_deg);
    const double margin   = margin_cfg > 0.0 ? margin_cfg : 0.5 * flat + MAGMA_IRON_RIM_CLEARANCE;
    CraterIron   iron;
    iron.start_R = crater_r + margin;
    iron.turns   = turns_cfg > 0
        ? turns_cfg
        : std::max(1, int(std::ceil(iron.start_R / std::max(0.1, MAGMA_IRON_STEP_FRACTION * flat))));
    return iron;
}

// Crater-iron speed: the explicit setting, or 0 = the region's ironing speed. Never travel
// speed -- plowing a molten rim at 150-300 mm/s tears it away instead of shearing it back in.
inline double resolve_iron_speed(double iron_speed_cfg, double ironing_speed)
{
    return iron_speed_cfg > 0.0 ? iron_speed_cfg : ironing_speed;
}

// Volumetric injection rate. 0 = the filament's max volumetric speed; an explicit rate is
// capped by it. Floored at 1 mm3/s: slower holds the nozzle in the cell long enough to melt it.
// The flags say which substitutions happened, so the emitter can report them.
struct InjectionRate
{
    double mm3_s    = 0.0;
    bool   invented = false;  // 0 was set but the filament declares no max volumetric speed
    bool   floored  = false;  // the resolved rate was below the floor
};
inline InjectionRate resolve_injection_rate(double speed_cfg, double max_vol)
{
    InjectionRate r;
    if (speed_cfg <= 0.0) {
        r.invented = max_vol <= 0.0;
        r.mm3_s    = r.invented ? 1.0 : max_vol;
    } else
        r.mm3_s = max_vol > 0.0 ? std::min(speed_cfg, max_vol) : speed_cfg;
    if (r.mm3_s < 1.0) {
        r.floored = true;
        r.mm3_s   = 1.0;
    }
    return r;
}

struct TriangleGeometry final : public MagmaGeometry
{
    // Insetting a 60deg corner by lw/2 shortens each side by lw*sqrt3.
    double open_edge_length(double spacing, double line_width) const override {
        return triangle_side_length(spacing) - line_width * SQRT3;
    }

    double inset_open_area(double spacing, double line_width) const override {
        return inset_triangle_area(spacing, line_width);
    }

    // Circumscribed circle of the inset triangle: 2*(side - lw*sqrt3)/sqrt3.
    double opening_diameter(double spacing, double line_width) const override {
        double inset_side = triangle_side_length(spacing) - line_width * SQRT3;
        return inset_side > 0.0 ? 2.0 * inset_side / SQRT3 : 0.0;
    }

    // With altitude spacing = interior + lw the inset side is (2*interior - lw)/sqrt3, so the
    // inradius is inset_side/(2*sqrt3) = (2*interior - lw)/6.
    double inscribed_radius(double interior_width, double line_width) const override {
        return std::max(0.0, (2.0 * interior_width - line_width) / 3.0) * 0.5;
    }

    // 3 line families crossing at 60deg per vertex: (3*sqrt3/4) * w^2.
    double vertex_overlap_excess_area(double line_width) const override {
        return (3.0 * SQRT3 / 4.0) * line_width * line_width;
    }

    double auto_window_height(double interior_width, double line_width) const override {
        const double spacing   = cell_spacing_from_geometry(interior_width, line_width);
        const double open_edge = open_edge_length(spacing, line_width);
        return open_edge > 0.0 ? inset_triangle_area(spacing, line_width) / open_edge : 0.0;
    }

    int max_neighbors() const override { return 3; }
};

inline const MagmaGeometry& triangle_geometry() {
    static const TriangleGeometry s_geom;
    return s_geom;
}

// Triangle cells are CellId (a, b, c): up triangles have a + b + c == 2, down triangles == 1.
// Lattice (lx, ly) maps to world as
//   px = lx * edge_length + ly * edge_length / 2,   py = ly * cell_spacing
// so rows are cell_spacing apart and each row shifts by half an edge.
class TriangleLattice : public MagmaLattice {
public:
    TriangleLattice() : m_cell_spacing(0), m_edge_length(0), m_offset_x(0), m_offset_y(0) {}

    explicit TriangleLattice(double cell_spacing, double offset_x = 0.0, double offset_y = 0.0)
        : m_cell_spacing(cell_spacing)
        , m_edge_length(triangle_side_length(cell_spacing))
        , m_offset_x(offset_x)
        , m_offset_y(offset_y)
    {}

    double cell_spacing() const override { return m_cell_spacing; }
    double edge_length() const override { return m_edge_length; }
    double offset_x() const override { return m_offset_x; }
    double offset_y() const override { return m_offset_y; }

    // ---- topology ----

    int max_neighbors() const override { return 3; }

    bool is_up(const CellId& cell) const override {
        return (cell.a + cell.b + cell.c) == 2;
    }

    // Up and down triangles neighbour each other: step one coordinate by -1 / +1.
    std::vector<CellId> neighbors(const CellId& cell) const override {
        const int a = cell.a, b = cell.b, c = cell.c;
        if (is_up(cell))
            return {{ {a-1,b,c}, {a,b-1,c}, {a,b,c-1} }};
        else
            return {{ {a+1,b,c}, {a,b+1,c}, {a,b,c+1} }};
    }

    // ---- coordinate transforms ----

    Vec2d to_world(double lx, double ly) const override {
        return Vec2d(
            lx * m_edge_length + ly * m_edge_length * 0.5 + m_offset_x,
            ly * m_cell_spacing + m_offset_y
        );
    }

    std::pair<double, double> to_lattice(double px, double py) const override {
        double adjusted_x = px - m_offset_x;
        double adjusted_y = py - m_offset_y;
        double ly = adjusted_y / m_cell_spacing;
        double lx = (adjusted_x - ly * m_edge_length * 0.5) / m_edge_length;
        return {lx, ly};
    }

    // ---- cell geometry ----

    CellId cell_at(double px, double py) const override {
        auto [lx, ly] = to_lattice(px, py);

        int col = static_cast<int>(std::floor(lx));
        int row = static_cast<int>(std::floor(ly));

        double fx = lx - col;
        double fy = ly - row;

        // Lower-left half of the lattice parallelogram is the up triangle.
        bool cell_is_up = (fx + fy) < 1.0;

        int c = cell_is_up ? (2 - col - row) : (1 - col - row);
        return CellId(col, row, c);
    }

    std::vector<Vec2d> cell_corners(const CellId& cell) const override;
    Vec2d cell_center(const CellId& cell) const override;
    std::vector<CellId> enumerate_cells(const BoundingBox& bbox) const override;

private:
    double m_cell_spacing;
    double m_edge_length;
    double m_offset_x;
    double m_offset_y;
};

// A window is the gap in the shared wall between the two cells of a U-tube.
struct WindowSpec {
    double window_height_mm = 0.4;

    WindowSpec() = default;

    static WindowSpec from_config(
        const MagmaGeometry& geom,
        float config_window_height_mm,      // 0 = auto
        float interior_width,
        float line_width,
        float layer_height                  // auto adds one layer
    );

};

} // namespace magma
} // namespace Slic3r

#endif // slic3r_Magma_MagmaTriangleCell_hpp_
