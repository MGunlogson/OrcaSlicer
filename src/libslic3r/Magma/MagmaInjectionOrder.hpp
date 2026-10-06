#ifndef slic3r_MagmaInjectionOrder_hpp_
#define slic3r_MagmaInjectionOrder_hpp_

#include "../Point.hpp"

#include <functional>
#include <vector>

namespace Slic3r {
namespace magma {

// Injection-phase timing, read once from the print config. Converts a visiting order into
// elapsed seconds, so heat decay accounts for travel time as cooling time.
struct InjectionTiming {
    double travel_speed_mm_s    = 150.0;  // XY travel speed between injections
    double per_injection_fixed_s = 0.0;   // ~constant per-injection time: z-hops + dwell (+ slam)
    double vol_speed_mm3_s      = 10.0;   // volumetric injection rate (extrude time = volume / rate)
};

// Visiting order over injection points: a permutation of indices into world_pts (scaled
// world-space XY).
//   spread_heat == false : travel-optimal order (chain_points).
//   spread_heat == true  : heat-spread order that keeps spatially-near injections far apart
//                          in time. Deterministic, O(n^2).
// injected_volume_mm3 is parallel to world_pts (fill factor applied) and sets each
// injection's duration.
std::vector<size_t> order_injection_points(
    const Points& world_pts,
    const std::vector<double>& injected_volume_mm3,
    const InjectionTiming& timing,
    bool spread_heat,
    const std::function<void()>& throw_if_canceled = {});

} // namespace magma
} // namespace Slic3r

#endif // slic3r_MagmaInjectionOrder_hpp_
