#include "MagmaSpiralOffset.hpp"

#include <cmath>
#include <algorithm>

namespace Slic3r {
namespace magma {

SpiralParams compute_spiral_params(float interior_width, float line_width, bool enabled)
{
    SpiralParams params;
    params.enabled = enabled;

    if (!enabled) {
        params.spiral_radius = 0.f;
        params.angle_per_layer = 0.f;
        return params;
    }

    // Constraint 1: Printability - 40% line overlap between layers
    constexpr float target_line_overlap = 0.40f;
    const float max_disp_line = (1.0f - target_line_overlap) * line_width;

    // Constraint 2: Tube continuity - 75% tube area overlap between layers
    constexpr float target_tube_overlap = 0.75f;
    const float max_disp_tube = (1.0f - target_tube_overlap) * interior_width;

    // Use the more restrictive constraint
    const float max_displacement = std::min(max_disp_line, max_disp_tube);

    // Interlock constraint: swept circles of adjacent tubes should touch
    const float cell_spacing = static_cast<float>(cell_spacing_from_geometry(interior_width, line_width));
    params.spiral_radius = cell_spacing / 2.0f;

    // Per-layer angle follows from radius and displacement constraints
    params.angle_per_layer = max_displacement / params.spiral_radius;

    return params;
}

Vec2d compute_spiral_offset(const SpiralParams &params, int layer_id)
{
    if (!params.enabled)
        return Vec2d(0.0, 0.0);

    const float layer_angle = float(layer_id) * params.angle_per_layer;
    return Vec2d(
        params.spiral_radius * std::cos(layer_angle),
        params.spiral_radius * std::sin(layer_angle)
    );
}

TriangleLattice lattice_for_layer(double cell_spacing, const SpiralParams &params, int layer_id)
{
    Vec2d offset = compute_spiral_offset(params, layer_id);
    return TriangleLattice(cell_spacing, offset.x(), offset.y());
}

} // namespace magma
} // namespace Slic3r
