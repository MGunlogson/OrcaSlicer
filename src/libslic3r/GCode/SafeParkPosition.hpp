#ifndef slic3r_GCode_SafeParkPosition_hpp_
#define slic3r_GCode_SafeParkPosition_hpp_

#include "../ExPolygon.hpp"
#include "../Point.hpp"

#include <optional>
#include <string>

namespace Slic3r {

class GCode;
class Layer;
class SupportLayer;

// Park position tier; lower = ooze is less likely to affect print quality.
enum class ParkPriority {
    Empty       = 1,  // Outside cumulative object footprint
    Support     = 2,  // Over support on current layer
    SparseInfill = 3, // Over sparse infill on current layer
    SolidInfill = 4,  // Over solid internal infill on current layer
    None        = 5   // No safe XY position found
};

struct ParkResult {
    std::optional<Point> position;                    // nullopt = z-hop only (Priority None)
    ParkPriority         priority = ParkPriority::None;

    bool needs_z_hop() const { return priority >= ParkPriority::SparseInfill; }
};

// E word for the park path's raw extra retract / unretract. Relative E emits the delta;
// absolute E needs a target, since "E-2" at E=850 would be an ~850mm retraction. Both are
// measured from e_retracted, the E the writer holds after the normal retract, and the
// unretract returns E there so the writer's own unretract still restores the rest.
inline double park_extra_retract_e(bool relative_e, double e_retracted, double extra)
{
    return relative_e ? -extra : e_retracted - extra;
}
inline double park_extra_unretract_e(bool relative_e, double e_retracted, double extra)
{
    return relative_e ? extra : e_retracted;
}

// Finds where to park the nozzle during a temperature change. Tracks the cumulative XY
// footprint of all objects printed so far, simplified to ~1mm since parking needs only
// coarse precision. Support and infill type are checked on the current layer only: support
// on lower layers is empty space at this height.
class SafeParkPosition {
public:
    // Add one object's layer slices to the cumulative footprint; call for each object at this Z.
    void update(const Layer* object_layer);

    // Safe XY park position near nozzle_pos, by priority:
    // empty > support > sparse infill > solid infill > none.
    ParkResult find_safe_position(
        const Layer* object_layer,
        const SupportLayer* support_layer,
        const Point& nozzle_pos,
        coord_t margin = scale_(2.0)) const;

    // Park and wait for the temperature:
    // retract -> (z-hop if needed) -> XY travel -> extra retract -> M109 -> extra unretract.
    static std::string park_and_set_temp(
        GCode& gcodegen,
        const ParkResult& park,
        double print_z,
        double park_z_hop,
        double extra_retract,
        int target_temp,
        const char* z_comment,
        const char* xy_comment);

private:
    ExPolygons m_cumulative_footprint;  // simplified cumulative object area

    // Find nearest centroid in safe regions to nozzle_pos.
    static std::optional<Point> nearest_centroid(
        const ExPolygons& safe_regions, const Point& nozzle_pos);

    static constexpr double SIMPLIFY_TOLERANCE = 1.0;  // mm
};

} // namespace Slic3r

#endif // slic3r_GCode_SafeParkPosition_hpp_
