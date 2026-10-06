#ifndef slic3r_FillMagma_hpp_
#define slic3r_FillMagma_hpp_

#include "FillBase.hpp"
#include "../Magma/MagmaTriangleCell.hpp"

namespace Slic3r {

namespace magma { class MagmaTubeMap; }

// Base for all Magma patterns. Fill.cpp attaches the object's tube map through this class,
// so subclasses only supply line generation.
class FillMagmaBase : public Fill
{
public:
    // Non-owning; set by Fill.cpp.
    const magma::MagmaTubeMap* tube_map = nullptr;

    // Throws if no tube map is attached, rather than printing the zone hollow.
    const magma::MagmaTubeMap& require_tube_map(const char *pattern_name) const;

    bool is_self_crossing() override { return true; }

    // Keep connect_infill's ordering; re-sorting scatters the many window-gap fragments.
    bool no_sort() const override { return true; }
};

// Triangle grid of three line families, translated per layer by the spiral offset. Window gaps
// are cut into the lines as they are generated.
class FillMagmaTriangle : public FillMagmaBase
{
public:
    Fill* clone() const override { return new FillMagmaTriangle(*this); }
    ~FillMagmaTriangle() override = default;

protected:
    // No rotation between layers; the spiral offset is a translation.
    float _layer_angle(size_t idx) const override { return 0.f; }

    // Grid is aligned to the world origin so cells line up across layers.
    std::pair<float, Point> _infill_direction(const Surface *surface) const override;

    void _fill_surface_single(
        const FillParams &params,
        unsigned int thickness_layers,
        const std::pair<float, Point> &direction,
        ExPolygon expolygon,
        Polylines &polylines_out) override;
};

// Square grid: one horizontal line per row and one vertical line per column, with window gaps.
class FillMagmaRectilinear : public FillMagmaBase
{
public:
    Fill* clone() const override { return new FillMagmaRectilinear(*this); }
    ~FillMagmaRectilinear() override = default;

protected:
    float _layer_angle(size_t idx) const override { return 0.f; }
    std::pair<float, Point> _infill_direction(const Surface *surface) const override;
    void _fill_surface_single(
        const FillParams &params,
        unsigned int thickness_layers,
        const std::pair<float, Point> &direction,
        ExPolygon expolygon,
        Polylines &polylines_out) override;
};

// Trihexagonal: hexagon hubs and triangle vents (MagmaTriHexCell.hpp), drawn as three line
// families offset half an index from the triangle pattern's. Windows are cut in the shared
// walls of open hub<->vent pairs.
class FillMagmaTriHex : public FillMagmaBase
{
public:
    Fill* clone() const override { return new FillMagmaTriHex(*this); }
    ~FillMagmaTriHex() override = default;

protected:
    float _layer_angle(size_t idx) const override { return 0.f; }
    std::pair<float, Point> _infill_direction(const Surface *surface) const override;
    void _fill_surface_single(
        const FillParams &params,
        unsigned int thickness_layers,
        const std::pair<float, Point> &direction,
        ExPolygon expolygon,
        Polylines &polylines_out) override;
};

// Pointy-top honeycomb (HexLattice), cells paired into 2-cell U-tubes. Drawn as continuous
// zigzags like Orca's honeycomb, so vertical walls are doubled and slants single; windows are
// cut by subtracting a rectangle over each open pair's shared wall.
class FillMagmaHoneycomb : public FillMagmaBase
{
public:
    Fill* clone() const override { return new FillMagmaHoneycomb(*this); }
    ~FillMagmaHoneycomb() override = default;

protected:
    float _layer_angle(size_t idx) const override { return 0.f; }
    std::pair<float, Point> _infill_direction(const Surface *surface) const override;
    void _fill_surface_single(
        const FillParams &params,
        unsigned int thickness_layers,
        const std::pair<float, Point> &direction,
        ExPolygon expolygon,
        Polylines &polylines_out) override;
};

} // namespace Slic3r

#endif // slic3r_FillMagma_hpp_
