#ifndef slic3r_FillMagma_hpp_
#define slic3r_FillMagma_hpp_

#include "FillRectilinear.hpp"

namespace Slic3r {

// Magma Triangle infill pattern for vertical reinforcement
//
// Creates a triangle grid pattern that shifts slightly each layer to create
// interlocking spirals. The shift is calculated to maintain ~30% overlap with
// the same-direction lines from 3 layers below (triangle has 3-layer cycle).
//
// Future enhancements:
// - Dynamic corner width for rounder tube cross-sections
// - Window gaps for U-tube pairing
// - Stagger levels for Z-offset weak plane avoidance
class FillMagmaTriangle : public FillRectilinear
{
public:
    Fill* clone() const override { return new FillMagmaTriangle(*this); }
    ~FillMagmaTriangle() override = default;

    Polylines fill_surface(const Surface *surface, const FillParams &params) override;

    // Triangle pattern is self-crossing (3 directions intersect)
    bool is_self_crossing() override { return true; }

protected:
    // Fixed angle - pattern doesn't rotate between layers
    // (spiral interlock comes from circular translation offset)
    float _layer_angle(size_t idx) const override { return 0.f; }

    // Calculate circular spiral offset for this layer
    // The pattern translates in a circle, creating helix-shaped tubes
    // that interlock when filled with plastic.
    // Rotation period is dynamically calculated for max spiral rate at 40% overlap.
    // Output: (offset_x, offset_y) in mm
    void calculate_spiral_offset(const FillParams &params, float interior_width, float &offset_x, float &offset_y) const;

    // Project (x,y) offset onto the perpendicular of a line direction
    // Returns the pattern_shift value for that line direction
    float project_offset_to_shift(float offset_x, float offset_y, float line_angle) const;
};

// Magma Hex infill pattern (future)
// Similar to triangle but uses hexagonal cells
class FillMagmaHex : public FillRectilinear
{
public:
    Fill* clone() const override { return new FillMagmaHex(*this); }
    ~FillMagmaHex() override = default;

    Polylines fill_surface(const Surface *surface, const FillParams &params) override;
    bool is_self_crossing() override { return true; }

protected:
    float _layer_angle(size_t idx) const override { return 0.f; }
};

} // namespace Slic3r

#endif // slic3r_FillMagma_hpp_
