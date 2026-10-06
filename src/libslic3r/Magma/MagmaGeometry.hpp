#ifndef slic3r_Magma_MagmaGeometry_hpp_
#define slic3r_Magma_MagmaGeometry_hpp_

#include <algorithm>

namespace Slic3r {
namespace magma {

// Minimum tube height: the shared wall above the window must be at least as tall as the window,
// or the pair fills as one cavity instead of a U. Shared by MagmaTubeMap::build (measured layer
// heights) and resolve_magma (nominal height, for warnings) so warnings match the g-code.
inline double magma_min_tube_height(double window_height_mm) {
    return 2.0 * window_height_mm;
}

// Fraction of the area-matched window height used by the auto window. The nozzle bore, not the
// window, is the flow restriction (a 1.6mm tube is ~9x a 0.6mm bore), so a half-size window still
// flows freely, and since min tube height is 2x the window it widens the legal tube-height band
// that stagger needs.
inline constexpr double MAGMA_AUTO_WINDOW_FRACTION = 0.5;

// Per-cell-shape geometry. Each Magma pattern implements this; the tube map, injection g-code
// and Print::validate() all consume it so every shape formula lives in one place.
//  - spacing:    center-to-center distance between adjacent parallel infill lines.
//  - line_width: deposited bead width, unless the caller deliberately passes a nominal one.
//  - Lengths in mm, areas in mm^2, volumes in mm^3; no scaled coordinates.
// Cell positions and neighbours belong to MagmaLattice; the nozzle seal math is shape-agnostic
// and only takes opening_diameter().
struct MagmaGeometry
{
    virtual ~MagmaGeometry() = default;

    // Width a window cut in a shared wall opens: the length of that edge left open after insetting
    // each wall by half the line width. Zero or less when the walls close the cell, which
    // Print::validate() rejects.
    virtual double open_edge_length(double spacing, double line_width) const = 0;

    // Open tube cross-section after insetting each wall by half the line width (mm^2).
    virtual double inset_open_area(double spacing, double line_width) const = 0;

    // Diameter of the circle that must be covered to seal the tube top
    // (circumscribed circle of the inset opening). Feeds the z-slam/seal math.
    virtual double opening_diameter(double spacing, double line_width) const = 0;

    // Radius of the largest circle inside the open tube (the usable bore). Takes line_width
    // because a triangle's inset shrinks the inradius by more than half a bead.
    virtual double inscribed_radius(double interior_width, double line_width) const = 0;

    // Excess area (mm^2) per cell where infill lines overlap at vertices; subtracted from the
    // injection volume.
    virtual double vertex_overlap_excess_area(double line_width) const = 0;

    // Window height (mm) at which the window's flow area (open edge x height) equals the tube's
    // open area. Zero when the cell has no open edge.
    virtual double auto_window_height(double interior_width, double line_width) const = 0;

    // Topology constants.
    virtual int max_neighbors() const = 0;   // edge-sharing neighbor count (3/4/6)
};

} // namespace magma
} // namespace Slic3r

#endif // slic3r_Magma_MagmaGeometry_hpp_
