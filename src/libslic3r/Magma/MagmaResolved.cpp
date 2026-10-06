#include "MagmaResolved.hpp"

#include <algorithm>

#include "../Flow.hpp"
#include "MagmaPatterns.hpp"
#include "MagmaTriangleCell.hpp"
#include "MagmaTubeMap.hpp"

namespace Slic3r {
namespace magma {


bool resolve_magma(const PrintRegionConfig &region,
                   const PrintObjectConfig &object,
                   const PrintConfig       &print,
                   MagmaResolved           &out)
{
    const InfillPattern pattern = magma_effective_pattern(region);
    if (! is_magma_pattern(pattern))
        return false;

    MagmaResolved r;
    r.pattern  = pattern;
    r.geometry = &magma_geometry_for(pattern);

    // Size the lattice against the extruder that prints it: with dual infill that is the outer
    // zone filament, since the lattice lives in the outer zone.
    // Ids are 1-based with 0 = "Default". PrintApply resolves 0 before slicing; the max() covers
    // the GUI readout, which resolves against an unresolved preset.
    const int lattice_filament = region.dual_infill_enabled.value
                                     ? region.dual_infill_outer_filament.value
                                     : region.sparse_infill_filament_id.value;
    r.sparse_extruder = std::max(0, lattice_filament - 1);
    r.nozzle_diameter = print.nozzle_diameter.get_at(r.sparse_extruder);

    r.line_width = region.sparse_infill_line_width.get_abs_value(r.nozzle_diameter);
    if (r.line_width <= 0.0)
        // 0 means "auto"; resolve it the way the slicer does, not as the nozzle diameter.
        r.line_width = Flow::auto_extrusion_width(frInfill, float(r.nozzle_diameter));

    // Seal settings are per-extruder and indexed by the injection extruder, the only tip that
    // touches a tube.
    const int inj_filament = object.magma_injection_filament.value;
    r.injection_extruder   = inj_filament > 0 ? inj_filament - 1 : r.sparse_extruder;
    r.injection_nozzle_diameter = print.nozzle_diameter.get_at(r.injection_extruder);
    r.injection_flat_is_estimate =
        object.magma_nozzle_outer_diameter.get_at(r.injection_extruder) <= 0.0;
    r.injection_nozzle_flat     = resolve_nozzle_flat(
        object.magma_nozzle_outer_diameter.get_at(r.injection_extruder), r.injection_nozzle_diameter);
    r.cone_half_angle_deg = object.magma_nozzle_cone_half_angle.get_at(r.injection_extruder);
    r.seal_press        = std::max(0.0, object.magma_seal_press.value);

    // The tube width is configured; the seal depth follows from it and the nozzle geometry.
    r.interior_width   = effective_interior_width(object.magma_interior_width.value);
    r.cell_spacing     = cell_spacing_from_geometry(r.interior_width, r.line_width);
    r.opening_diameter = r.geometry->opening_diameter(r.cell_spacing, r.line_width);
    r.bore_diameter    = 2.0 * r.geometry->inscribed_radius(r.interior_width, r.line_width);
    r.seal_depth       = auto_seal_depth(r.opening_diameter, r.injection_nozzle_flat,
                                         r.cone_half_angle_deg, r.seal_press);

    // Tube height bounds, through the same helpers MagmaTubeMap::build() uses.
    {
        const double nominal_lh = std::max(0.0, object.layer_height.value);
        const WindowSpec ws = WindowSpec::from_config(
            *r.geometry,
            static_cast<float>(object.magma_window_height_mm.value),
            static_cast<float>(r.interior_width),
            static_cast<float>(r.line_width),
            static_cast<float>(nominal_lh));
        r.min_tube_height = magma_min_tube_height(ws.window_height_mm);
        r.max_tube_height = std::max(std::max(1.0, object.magma_tube_height.value),
                                     r.min_tube_height);
    }

    // Plunge is added on top of the seal depth, bounded only by clamp_plunge_depth().
    r.plunge_requested = requested_plunge_depth(object.magma_injection_plunge.value,
                                                object.magma_injection_plunge_depth.value);
    r.plunge_depth     = clamp_plunge_depth(r.seal_depth, r.plunge_requested);

    r.grip             = corner_grip(r.seal_press, r.plunge_depth, r.cone_half_angle_deg);
    r.cone_at_seal     = cone_diameter_at(r.seal_depth, r.injection_nozzle_flat, r.cone_half_angle_deg);

    out = r;
    return true;
}

} // namespace magma
} // namespace Slic3r
