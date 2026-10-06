#include <catch2/catch_all.hpp>

#include <cmath>
#include <utility>
#include <vector>

#include "libslic3r/libslic3r.h"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/GCode/SafeParkPosition.hpp"
#include "libslic3r/Magma/MagmaGeometry.hpp"
#include "libslic3r/Magma/MagmaPatterns.hpp"
#include "libslic3r/Magma/MagmaResolved.hpp"
#include "libslic3r/Magma/MagmaTriangleCell.hpp"
#include "libslic3r/Magma/MagmaTubeMap.hpp"
#include "libslic3r/Magma/MagmaTubeSolver.hpp"

using namespace Slic3r;
using namespace Slic3r::magma;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

struct PatternCase {
    const char   *name;
    InfillPattern pattern;
    // Bore (inscribed circle) over seal opening (circumscribed circle): cos(pi/n) for a regular n-gon.
    double        bore_over_opening;
};

const std::vector<PatternCase> &all_patterns()
{
    static const std::vector<PatternCase> cases = {
        { "MagmaTriangle",    ipMagmaTriangle,    0.5                  },
        { "MagmaRectilinear", ipMagmaRectilinear, 1.0 / std::sqrt(2.0) },
        { "MagmaHoneycomb",   ipMagmaHoneycomb,   std::sqrt(3.0) / 2.0 },
        { "MagmaTriHex",      ipMagmaTriHex,      std::sqrt(3.0) / 2.0 },
    };
    return cases;
}

// Representative geometries: a small tube on a fine nozzle, a typical one, a large one.
struct Sizing { double interior; double line_width; };
const std::vector<Sizing> &sizings()
{
    static const std::vector<Sizing> s = {
        { 1.20, 0.45 }, { 1.62, 0.60 }, { 2.50, 0.66 }, { 4.00, 0.80 },
    };
    return s;
}

// A classic E3D 0.6mm nozzle: 1.70mm flat, 30 degree cone.
constexpr double REF_FLAT = 1.70, REF_CONE = 30.0, REF_NOZZLE = 0.60;
constexpr double REF_LINE_WIDTH = 0.42, REF_INTERIOR = 1.6;
constexpr double REF_PRESS = 0.05, REF_PLUNGE = 0.05;

// Resolves through resolve_magma(), the entry point MagmaTubeMap::build, Print::validate and
// PresetHints use, so the test covers the shipped chain rather than a reimplementation.
MagmaResolved resolve_ref(InfillPattern pattern)
{
    PrintRegionConfig region;
    PrintObjectConfig object;
    PrintConfig       printer;

    printer.nozzle_diameter.values = { REF_NOZZLE };

    region.sparse_infill_pattern.value      = pattern;
    region.sparse_infill_filament_id.value  = 1;
    region.dual_infill_enabled.value        = false;
    region.sparse_infill_line_width.value   = REF_LINE_WIDTH;
    region.sparse_infill_line_width.percent = false;

    object.magma_interior_width.value          = REF_INTERIOR;
    object.magma_nozzle_outer_diameter.values  = { REF_FLAT };
    object.magma_nozzle_cone_half_angle.values = { REF_CONE };
    object.magma_seal_press.value              = REF_PRESS;
    object.magma_injection_plunge.value        = true;
    object.magma_injection_plunge_depth.value  = REF_PLUNGE;

    MagmaResolved m;
    REQUIRE(resolve_magma(region, object, printer, m));
    return m;
}

} // namespace

TEST_CASE("The bore is the cell's inscribed circle, not the interior width", "[Magma][Geometry]")
{
    // interior_width * 0.5 is right for square, hex and tri-hex but not triangle.
    for (const PatternCase &pc : all_patterns()) {
        const MagmaGeometry &geom = magma_geometry_for(pc.pattern);
        for (const Sizing &sz : sizings()) {
            const double spacing = cell_spacing_from_geometry(sz.interior, sz.line_width);
            const double opening = geom.opening_diameter(spacing, sz.line_width);
            const double bore    = 2.0 * geom.inscribed_radius(sz.interior, sz.line_width);
            INFO(pc.name << " interior=" << sz.interior << " lw=" << sz.line_width
                         << " opening=" << opening << " bore=" << bore);
            CHECK(bore / opening == Catch::Approx(pc.bore_over_opening).epsilon(1e-6));
        }
    }
}

TEST_CASE("The seal press is depth past first contact, whatever the tube", "[Magma][Seal]")
{
    // Includes openings the flat already spans, where first contact is at the surface.
    const double cone = 30.0;
    for (double press : { 0.0, 0.05, 0.0866, 0.3 })
        for (double flat : { 1.2, 1.7, 2.4 })
            for (double opening : { 1.0, 2.0, 2.9, 3.6 }) {
                const double d_seal    = auto_seal_depth(opening, flat, cone, press);
                const double d_contact = seal_depth_for_opening(opening, flat, cone);
                INFO("flat=" << flat << " opening=" << opening << " press=" << press);
                CHECK(d_seal - d_contact == Catch::Approx(press).epsilon(1e-9));
            }
}

TEST_CASE("The effective pattern follows the dual-infill state", "[Magma][Geometry]")
{
    PrintRegionConfig region;
    region.sparse_infill_pattern.value     = ipMagmaRectilinear;
    region.dual_infill_outer_pattern.value = ipMagmaHoneycomb;

    region.dual_infill_enabled.value = false;
    CHECK(magma_effective_pattern(region) == ipMagmaRectilinear);   // dual off: sparse wins
    region.dual_infill_enabled.value = true;
    CHECK(magma_effective_pattern(region) == ipMagmaHoneycomb);     // dual on: outer wins
    // A non-Magma outer pattern (a zone with no injectable channels) is returned as chosen.
    region.dual_infill_outer_pattern.value = ipGrid;
    CHECK(magma_effective_pattern(region) == ipGrid);
}

TEST_CASE("The auto window matches the tube's open area and fits its open edge", "[Magma][Geometry]")
{
    const double layer_height = 0.2;
    for (const PatternCase &pc : all_patterns()) {
        const MagmaGeometry &geom = magma_geometry_for(pc.pattern);
        for (const Sizing &sz : sizings()) {
            const double spacing   = cell_spacing_from_geometry(sz.interior, sz.line_width);
            const double open_edge = geom.open_edge_length(spacing, sz.line_width);
            INFO(pc.name << " interior=" << sz.interior << " lw=" << sz.line_width);
            REQUIRE(open_edge > 0.0);
            // Tri-hex matches the vent and floors the result, so only the others are area / edge.
            if (pc.pattern != ipMagmaTriHex)
                CHECK_THAT(geom.auto_window_height(sz.interior, sz.line_width),
                           WithinRel(geom.inset_open_area(spacing, sz.line_width) / open_edge, 1e-9));
            // The auto window is never taller than the edge it is cut in, plus the one added layer.
            const WindowSpec ws = WindowSpec::from_config(geom, 0.f, float(sz.interior), float(sz.line_width),
                                                          float(layer_height));
            CHECK(ws.window_height_mm <= open_edge + layer_height + 1e-6);
        }
    }
}

TEST_CASE("Honeycomb geometry and lattice agree on the open edge", "[Magma][Geometry]")
{
    const MagmaGeometry &geom = magma_geometry_for(ipMagmaHoneycomb);
    for (const Sizing &sz : sizings()) {
        const double spacing = cell_spacing_from_geometry(sz.interior, sz.line_width);
        const auto   lattice = make_magma_lattice(ipMagmaHoneycomb, spacing, 0.0, 0.0, sz.line_width);
        INFO("interior=" << sz.interior << " lw=" << sz.line_width);
        CHECK_THAT(lattice->edge_length(), WithinRel(geom.open_edge_length(spacing, sz.line_width), 1e-9));
    }
}

TEST_CASE("A cell whose walls meet has no open edge and no auto window", "[Magma][Geometry]")
{
    // A 0.2 mm tube in 0.6 mm lines: every shape's walls close it.
    const double interior = 0.2, lw = 0.6;
    for (const PatternCase &pc : all_patterns()) {
        const MagmaGeometry &geom = magma_geometry_for(pc.pattern);
        INFO(pc.name);
        if (geom.open_edge_length(cell_spacing_from_geometry(interior, lw), lw) <= 0.0)
            CHECK_THAT(geom.auto_window_height(interior, lw), WithinAbs(0.0, 1e-12));
    }
    CHECK(magma_geometry_for(ipMagmaTriangle).open_edge_length(cell_spacing_from_geometry(interior, lw), lw) <= 0.0);
}

TEST_CASE("Seal depth is solved for the cell's furthest corner", "[Magma][Seal]")
{
    // The cone must cover the circumscribed circle (the corners are reached last), then descend by the press.
    const double flat = 1.70, cone = 30.0, lw = 0.60;
    for (const PatternCase &pc : all_patterns()) {
        const MagmaGeometry &geom = magma_geometry_for(pc.pattern);
        for (double interior : { 1.2, 1.6, 2.0 }) {
            for (double press : { 0.0, 0.05, 0.2 }) {
                const double spacing = cell_spacing_from_geometry(interior, lw);
                const double opening = geom.opening_diameter(spacing, lw);
                const double d       = auto_seal_depth(opening, flat, cone, press);
                INFO(pc.name << " interior=" << interior << " press=" << press);
                CHECK(cone_diameter_at(d, flat, cone) + 1e-9 >= opening);
                CHECK(d - seal_depth_for_opening(opening, flat, cone)
                      == Catch::Approx(press).epsilon(1e-9));
                CHECK(2.0 * geom.inscribed_radius(interior, lw)
                      == Catch::Approx(opening * pc.bore_over_opening).epsilon(1e-9));
            }
        }
    }
}

TEST_CASE("Corner grip depends on press and plunge, nothing else", "[Magma][Seal]")
{
    // Grip is independent of tube size and seal depth; press and plunge are both depth past first
    // contact, so they contribute identically.
    const double cone = 30.0, T = std::tan(cone * MAGMA_DEG2RAD);
    for (double press : { 0.0, 0.05, 0.2 })
        for (double pl : { 0.0, 0.1, 0.3 }) {
            INFO("press=" << press << " plunge=" << pl);
            CHECK(corner_grip(press, pl, cone) == Catch::Approx((press + pl) * T).epsilon(1e-9));
        }
    // Covered but not gripped.
    CHECK(corner_grip(0.0, 0.0, cone) == Catch::Approx(0.0).margin(1e-12));
}

TEST_CASE("The plunge adds to the seal depth and only the slam clamp bounds the total", "[Magma][Seal]")
{
    for (double interior : { 1.2, 1.6, 2.0 })
        CHECK(effective_interior_width(interior) == Catch::Approx(interior));

    for (double seal : { 0.2, 0.6, 1.0 })
        for (double pl : { 0.0, 0.05, 0.3, 2.0 }) {
            const double got = clamp_plunge_depth(seal, pl);
            INFO("seal=" << seal << " plunge=" << pl);
            CHECK(got <= pl + 1e-9);
            CHECK(seal + got <= MAGMA_SLAM_CLAMP + 1e-9);
        }

    CHECK_THAT(requested_plunge_depth(false, 0.4), WithinAbs(0.0, 1e-12));
    CHECK_THAT(requested_plunge_depth(true, 0.4), WithinAbs(0.4, 1e-12));
    CHECK_THAT(requested_plunge_depth(true, -1.0), WithinAbs(0.0, 1e-12));
}

TEST_CASE("The resolver derives one consistent geometry per pattern", "[Magma][Seal]")
{
    for (const PatternCase &pc : all_patterns()) {
        const MagmaResolved m = resolve_ref(pc.pattern);
        INFO(pc.name);
        CHECK(m.line_width     == Catch::Approx(REF_LINE_WIDTH));
        CHECK(m.interior_width == Catch::Approx(REF_INTERIOR));
        CHECK(m.bore_diameter
              == Catch::Approx(m.opening_diameter * pc.bore_over_opening).epsilon(1e-9));
        // At seal depth the cone has passed first contact with the opening by `press` (a depth),
        // so it overlaps the opening by 2 * press * tan(theta).
        CHECK(cone_diameter_at(m.seal_depth, m.injection_nozzle_flat, m.cone_half_angle_deg)
              == Catch::Approx(m.opening_diameter
                               + 2.0 * m.seal_press
                                     * std::tan(m.cone_half_angle_deg * MAGMA_DEG2RAD))
                     .epsilon(1e-9));
        CHECK(m.total_depth() == Catch::Approx(m.seal_depth + m.plunge_depth).epsilon(1e-9));
        CHECK_FALSE(m.plunge_clamped());
    }
}

TEST_CASE("The crater iron starts with the nozzle flat outside the crater", "[Magma][Injection]")
{
    for (double flat : { 1.2, 1.7, 2.4 })
        for (double depth : { 0.0, 0.4, 1.0 }) {
            const CraterIron iron = plan_crater_iron(flat, REF_CONE, depth, /*margin*/ 0.0, /*turns*/ 0);
            const double crater_r = 0.5 * cone_diameter_at(depth, flat, REF_CONE);
            INFO("flat=" << flat << " depth=" << depth);
            // Auto margin: the flat's inner edge clears the crater rim by the rim clearance.
            CHECK_THAT(iron.start_R - 0.5 * flat, WithinAbs(crater_r + MAGMA_IRON_RIM_CLEARANCE, 1e-9));
            // Auto turns step inward by MAGMA_IRON_STEP_FRACTION of the flat per turn.
            CHECK(iron.turns == int(std::ceil(iron.start_R / (MAGMA_IRON_STEP_FRACTION * flat))));
            CHECK(iron.length(flat) > 2.0 * PI * iron.start_R);
        }

    // Explicit settings are used as given.
    const CraterIron set = plan_crater_iron(1.7, REF_CONE, 0.5, /*margin*/ 1.0, /*turns*/ 2);
    CHECK_THAT(set.start_R, WithinAbs(0.5 * cone_diameter_at(0.5, 1.7, REF_CONE) + 1.0, 1e-9));
    CHECK(set.turns == 2);
}

TEST_CASE("The injection rate resolves the way the settings describe", "[Magma][Injection]")
{
    // 0 = the filament's max volumetric speed.
    const InjectionRate from_filament = resolve_injection_rate(0.0, 12.0);
    CHECK_THAT(from_filament.mm3_s, WithinAbs(12.0, 1e-12));
    CHECK_FALSE(from_filament.invented);
    CHECK_FALSE(from_filament.floored);

    // An explicit rate is capped by the filament's max.
    CHECK_THAT(resolve_injection_rate(20.0, 12.0).mm3_s, WithinAbs(12.0, 1e-12));
    CHECK_THAT(resolve_injection_rate(8.0, 12.0).mm3_s, WithinAbs(8.0, 1e-12));
    CHECK_THAT(resolve_injection_rate(8.0, 0.0).mm3_s, WithinAbs(8.0, 1e-12));

    // 0 with no filament max is a substitution the emitter has to report.
    const InjectionRate invented = resolve_injection_rate(0.0, 0.0);
    CHECK(invented.invented);
    CHECK_THAT(invented.mm3_s, WithinAbs(1.0, 1e-12));

    // Below 1 mm3/s is floored, and says so.
    const InjectionRate floored = resolve_injection_rate(0.5, 12.0);
    CHECK(floored.floored);
    CHECK_THAT(floored.mm3_s, WithinAbs(1.0, 1e-12));

    // Iron speed: 0 borrows the ironing speed.
    CHECK_THAT(resolve_iron_speed(0.0, 25.0), WithinAbs(25.0, 1e-12));
    CHECK_THAT(resolve_iron_speed(40.0, 25.0), WithinAbs(40.0, 1e-12));
}

TEST_CASE("A print_z lookup resolves to the closest layer within EPSILON", "[Magma][Layers]")
{
    const std::vector<double> zs = { 0.2, 0.4, 0.40004, 0.6 };
    auto lookup = [&](double z) { return closest_print_z(zs.begin(), zs.end(), z, [](double v) { return v; }); };

    const auto [exact, d_exact] = lookup(0.6);
    REQUIRE(exact != zs.end());
    CHECK_THAT(*exact, WithinAbs(0.6, 1e-12));
    CHECK_THAT(d_exact, WithinAbs(0.0, 1e-12));

    // Layers merged within EPSILON: the nearer one wins.
    const auto [merged, d_merged] = lookup(0.40003);
    REQUIRE(merged != zs.end());
    CHECK_THAT(*merged, WithinAbs(0.40004, 1e-12));
    CHECK(d_merged < EPSILON);

    // Between layers: the closest is returned with its distance, for the caller to reject.
    const auto [between, d_between] = lookup(0.5);
    REQUIRE(between != zs.end());
    CHECK(d_between > EPSILON);

    // Above every layer.
    CHECK(lookup(1.0).first == zs.end());
}

TEST_CASE("The park retract is measured from the retracted E", "[Magma][Regression]")
{
    // A relative delta under absolute E would retract to E=-2 (~850 mm here). Under absolute E both
    // moves are measured from the already-retracted E, so the extra retract goes further back and
    // the unretract stops where the normal retract left it. Shared with ooze prevention.
    const double e_retracted = 849.2, extra = 2.0;

    CHECK_THAT(park_extra_retract_e  (true,  e_retracted, extra), WithinAbs(-2.0, 1e-9));
    CHECK_THAT(park_extra_unretract_e(true,  e_retracted, extra), WithinAbs( 2.0, 1e-9));
    CHECK_THAT(park_extra_retract_e  (false, e_retracted, extra), WithinAbs(847.2, 1e-9));
    CHECK_THAT(park_extra_unretract_e(false, e_retracted, extra), WithinAbs(849.2, 1e-9));
}

TEST_CASE("A raft's placeholder layers never become the minimum layer height", "[Magma][Regression]")
{
    // Layer data is indexed by absolute Layer::id(), so a raft leaves zero-height rows below the first
    // object layer. A zero minimum would make max_tube_height / min_height infinite.
    std::vector<LayerData> layers(8);
    for (auto &ld : layers) ld.height = 0.0;   // raft placeholders
    layers[3].height = 0.20;
    layers[4].height = 0.20;
    layers[5].height = 0.12;                   // the real minimum
    layers[6].height = 0.20;
    layers[7].height = 0.20;

    CHECK(min_positive_layer_height(layers, 3) == Catch::Approx(0.12));
    // Zeros are missing data, not thin layers.
    CHECK(min_positive_layer_height(layers, 0) == Catch::Approx(0.12));

    // No positive height: 0, so the caller can bail.
    std::vector<LayerData> all_zero(4);
    for (auto &ld : all_zero) ld.height = 0.0;
    CHECK(min_positive_layer_height(all_zero, 0) == Catch::Approx(0.0));
}

TEST_CASE("Honeycomb has no vertex overlap and tri-hex does", "[Magma][Regression]")
{
    // Honeycomb vertices are degree-3 line ends with no crossing; tri-hex lines cross. A spurious
    // overlap would thin the walls and open a leak path for the injected melt.
    const double lw = 0.42;
    CHECK(magma_geometry_for(ipMagmaHoneycomb).vertex_overlap_excess_area(lw)
          == Catch::Approx(0.0).margin(1e-12));
    CHECK(magma_geometry_for(ipMagmaTriHex).vertex_overlap_excess_area(lw) > 0.0);
}
