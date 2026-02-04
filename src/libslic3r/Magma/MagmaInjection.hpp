#ifndef slic3r_MagmaInjection_hpp_
#define slic3r_MagmaInjection_hpp_

#include "../Point.hpp"

#include <string>
#include <vector>

namespace Slic3r {

class GCode;
class Print;

namespace magma {

class MagmaTubeMap;

struct InjectionPoint {
    Vec2d  position;      // XY center of injection cell (mm, unscaled)
    double volume_mm3;    // volume to inject (after fill_factor)
    int    pair_index;    // index into tube map pairs
};

// Collect injection points for tubes whose cap layer == layer_id.
std::vector<InjectionPoint> collect_injection_points(
    const MagmaTubeMap& tube_map,
    int layer_id);

// Generate injection G-code for this layer.
// Returns empty string if no injection points exist.
std::string generate_injection_gcode(
    GCode& gcodegen,
    const MagmaTubeMap& tube_map,
    const std::vector<InjectionPoint>& points,
    double layer_z);

} // namespace magma
} // namespace Slic3r

#endif // slic3r_MagmaInjection_hpp_
