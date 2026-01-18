#include "FillMagma.hpp"

#include <cmath>
#include <boost/log/trivial.hpp>

namespace Slic3r {

// Calculate circular spiral offset for interlocking triangles
//
// The pattern shifts in a CIRCLE each layer, creating actual spirals.
// Adjacent tubes spiral around each other, mechanically interlocking when filled.
//
// Key insight: We want to spiral as FAST as possible while maintaining 40% overlap.
//
// Approach:
// 1. Radius is constrained by geometry (can't exceed hole or collide with neighbors)
// 2. Given fixed radius, back-calculate the maximum angle change that maintains 40% overlap
// 3. This gives the minimum rotation period (fastest possible spiral)
//
void FillMagmaTriangle::calculate_spiral_offset(const FillParams &params,
                                                 float interior_width,
                                                 float &offset_x, float &offset_y) const
{
    const float line_width = params.flow.width();

    // Spiral parameters
    constexpr float target_overlap = 0.40f;      // 40% overlap (overhang threshold)
    constexpr int cycle_layers = 3;              // Triangle repeats every 3 layers

    // 1. Calculate spiral radius from geometry constraints
    // Can't exceed half the hole, can't collide with neighbors
    const float cell_spacing = interior_width + line_width;
    const float spiral_radius = std::min(
        interior_width * 0.5f,                   // Don't exceed half the hole
        (cell_spacing - line_width) * 0.5f      // Don't collide with neighbors
    );

    // 2. Calculate max distance between same-direction lines for 40% overlap
    const float target_distance = (1.0f - target_overlap) * line_width;  // 0.6 * line_width

    // 3. Back-calculate angle change from radius and target distance
    // For circular motion: distance = 2 * radius * sin(angle_change / 2)
    // Solving: sin(angle_change / 2) = target_distance / (2 * radius)
    // But clamp to valid range for arcsin (can't exceed 1.0)
    const float sin_half_angle = std::min(1.0f, target_distance / (2.0f * spiral_radius));
    const float angle_change = 2.0f * std::asin(sin_half_angle);

    // 4. Calculate rotation period from angle change
    // rotation_period = cycle_layers * 2π / angle_change
    // (This is the minimum rotation period that maintains overlap)
    const float rotation_period = float(cycle_layers) * 2.0f * float(M_PI) / angle_change;

    // 5. Calculate current layer's position on the spiral circle
    const float layer_angle = float(this->layer_id) * 2.0f * float(M_PI) / rotation_period;

    offset_x = spiral_radius * std::cos(layer_angle);
    offset_y = spiral_radius * std::sin(layer_angle);
}

// Convert (x, y) offset to pattern_shift for a line at given angle
// pattern_shift moves lines perpendicular to their direction
float FillMagmaTriangle::project_offset_to_shift(float offset_x, float offset_y, float line_angle) const
{
    // Perpendicular direction to line_angle is (line_angle + 90°)
    // pattern_shift = dot((offset_x, offset_y), perpendicular_unit_vector)
    // perpendicular = (-sin(line_angle), cos(line_angle))
    return -offset_x * std::sin(line_angle) + offset_y * std::cos(line_angle);
}

Polylines FillMagmaTriangle::fill_surface(const Surface *surface, const FillParams &params)
{
    Polylines polylines_out;

    // Magma interior width: max dimension of triangular hole
    // For injection to work, nozzle must seal against the tube opening.
    // interior_width = nozzle_diameter + 0.2mm for sealing
    // TODO: Get from config, for now hardcode based on typical 0.4mm nozzle
    constexpr float interior_width = 0.6f;  // mm (nozzle + 0.2)

    const float line_width = params.flow.width();  // ~0.45mm

    // Line spacing = interior_width + line_width
    // This gives the center-to-center distance between parallel lines
    const float target_line_spacing = interior_width + line_width;

    // Save and override this->spacing
    // fill_surface_by_multilines uses formula: line_spacing = spacing * multiline / density
    // After density /= num_directions: effective line_spacing = spacing * multiline * num_directions / density
    // With density=1.0, multiline=1, num_directions=3: line_spacing = spacing * 3
    // So: spacing = target_line_spacing / 3
    double original_spacing = this->spacing;
    const_cast<FillMagmaTriangle*>(this)->spacing = target_line_spacing / 3.0;

    // Use full density (100%) - Magma doesn't use density concept
    FillParams magma_params = params;
    magma_params.density = 1.0f;

    // Calculate circular spiral offset for this layer
    float offset_x, offset_y;
    this->calculate_spiral_offset(magma_params, interior_width, offset_x, offset_y);

    // Triangle pattern: 3 sweeps at 0°, 60°, 120°
    // Each direction gets a different pattern_shift based on projecting
    // the circular (x,y) offset onto that direction's perpendicular axis.
    // This creates true circular motion of the triangle intersection points.
    constexpr float angle_0   = 0.f;
    constexpr float angle_60  = float(M_PI / 3.);
    constexpr float angle_120 = float(2. * M_PI / 3.);

    const float shift_0   = this->project_offset_to_shift(offset_x, offset_y, angle_0);
    const float shift_60  = this->project_offset_to_shift(offset_x, offset_y, angle_60);
    const float shift_120 = this->project_offset_to_shift(offset_x, offset_y, angle_120);

    if (magma_params.multiline > 1) {
        // Thick infill mode - use trapezoidal geometry
        // For trapezoidal, we use 0° and 90° directions
        const float shift_90 = this->project_offset_to_shift(offset_x, offset_y, float(M_PI / 2.));
        if (!this->fill_surface_trapezoidal(
                surface, magma_params,
                { { 0.f, shift_0 }, { float(M_PI / 2.), shift_90 } },
                polylines_out, 1))
            BOOST_LOG_TRIVIAL(error) << "FillMagmaTriangle::fill_surface_trapezoidal() failed.";
    } else {
        // Standard triangle: 3 line directions at 60° intervals
        // Each direction gets its projected shift for circular motion
        if (!this->fill_surface_by_multilines(
                surface, magma_params,
                { { angle_0, shift_0 },
                  { angle_60, shift_60 },
                  { angle_120, shift_120 } },
                polylines_out))
            BOOST_LOG_TRIVIAL(error) << "FillMagmaTriangle::fill_surface() failed to fill a region.";
    }

    // Restore original spacing
    const_cast<FillMagmaTriangle*>(this)->spacing = original_spacing;

    return polylines_out;
}

// FillMagmaHex - placeholder implementation
// Currently delegates to standard honeycomb, will be enhanced with Magma features
Polylines FillMagmaHex::fill_surface(const Surface *surface, const FillParams &params)
{
    Polylines polylines_out;

    // TODO: Implement hex-specific Magma pattern
    // For now, use standard triangle as placeholder
    // (Hex requires different geometry - will implement separately)

    BOOST_LOG_TRIVIAL(warning) << "FillMagmaHex: Using triangle pattern as placeholder";

    if (!this->fill_surface_by_multilines(
            surface, params,
            { { 0.f, 0.f }, { float(M_PI / 3.), 0.f }, { float(2. * M_PI / 3.), 0.f } },
            polylines_out))
        BOOST_LOG_TRIVIAL(error) << "FillMagmaHex::fill_surface() failed to fill a region.";

    return polylines_out;
}

} // namespace Slic3r
