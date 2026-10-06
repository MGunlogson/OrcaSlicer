#ifndef slic3r_MagmaInjection_hpp_
#define slic3r_MagmaInjection_hpp_

#include "../Point.hpp"

#include <string>
#include <vector>

namespace Slic3r {

class GCode;

namespace magma {

class MagmaTubeMap;

struct InjectionPoint {
    Vec2d  position;      // XY center of injection cell (mm, unscaled)
    double volume_mm3;    // measured tube volume; fill_factor is applied at injection
    int    pair_index;    // index into tube map pairs
    int    window_center_layer;  // center layer of window gap (for visualization)
};

// Ways an injection comes out weaker than the preview shows. Each leaves part of the lattice
// unfilled without any visible failure, so they are counted here for the export to report
// to the user rather than only logged.
struct InjectionDiagnostics {
    int tubes_no_volume        = 0;  // measured cavity was empty; that U-tube prints hollow
    int layers_no_fill_factor  = 0;  // fill factor <= 0; every injection on the layer dropped
    int injections_invented_speed = 0;  // no volumetric speed to inject at; a rate nobody chose
    int injections_speed_floored  = 0;  // asked for < 1 mm3/s; raised to the floor
    int injections_unlabelled  = 0;  // outside any exclude-object block; runs even if cancelled

    bool any() const
    {
        return tubes_no_volume || layers_no_fill_factor
            || injections_invented_speed || injections_speed_floored
            || injections_unlabelled;
    }
    void merge(const InjectionDiagnostics &o)
    {
        tubes_no_volume            += o.tubes_no_volume;
        layers_no_fill_factor      += o.layers_no_fill_factor;
        injections_invented_speed  += o.injections_invented_speed;
        injections_speed_floored   += o.injections_speed_floored;
        injections_unlabelled      += o.injections_unlabelled;
    }
    // One line per problem that actually occurred; empty when none did. Localized.
    std::vector<std::string> messages() const;
};

// Generate injection + crater-iron G-code for one object's points. Temperature and fan
// are handled by the caller (GCode.cpp injection phase). `diag` accumulates across calls.
std::string generate_injection_gcode(
    GCode& gcodegen,
    const MagmaTubeMap& tube_map,
    const std::vector<InjectionPoint>& points,
    double print_z,
    InjectionDiagnostics& diag);

// Park at a safe position (if enabled) and change nozzle temperature, waiting until reached.
// Used by GCode.cpp for the injection phase heat-up / cool-down, once per layer.
std::string park_and_set_temp(
    GCode& gcodegen,
    bool park_enabled,
    double print_z,
    double park_z_hop,
    double extra_retract,
    int target_temp,
    const char* z_comment,
    const char* xy_comment);

} // namespace magma
} // namespace Slic3r

#endif // slic3r_MagmaInjection_hpp_
