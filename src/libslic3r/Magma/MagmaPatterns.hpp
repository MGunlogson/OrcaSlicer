#ifndef slic3r_Magma_MagmaPatterns_hpp_
#define slic3r_Magma_MagmaPatterns_hpp_

#include "../PrintConfig.hpp"        // InfillPattern
#include "MagmaGeometry.hpp"
#include "MagmaLattice.hpp"
#include "MagmaTriangleCell.hpp"     // triangle_geometry(), TriangleLattice
#include "MagmaRectilinearCell.hpp"  // square_geometry(), RectilinearLattice
#include "MagmaTriHexCell.hpp"       // trihex_geometry(), TriHexLattice
#include "MagmaHexCell.hpp"          // hexagon_geometry(), HexLattice

#include <memory>

namespace Slic3r {
namespace magma {

// Maps an InfillPattern to its MagmaGeometry (scalar shape formulas) and MagmaLattice
// (topology and coordinates). A new Magma pattern adds a case to both selectors.

// Shape formulas for a pattern (shared, stateless instance).
inline const MagmaGeometry& magma_geometry_for(InfillPattern p)
{
    switch (p) {
    case ipMagmaRectilinear:
        return square_geometry();
    case ipMagmaTriHex:
        return trihex_geometry();
    case ipMagmaHoneycomb:
        return hexagon_geometry();
    case ipMagmaTriangle:
    default:
        return triangle_geometry();
    }
}

// Lattice for a pattern with the given cell spacing and spiral offset. line_width is used
// only by the hexagon lattice.
inline std::unique_ptr<MagmaLattice> make_magma_lattice(
    InfillPattern p, double cell_spacing, double offset_x, double offset_y, double line_width = 0.0)
{
    switch (p) {
    case ipMagmaRectilinear:
        return std::make_unique<RectilinearLattice>(cell_spacing, offset_x, offset_y);
    case ipMagmaTriHex:
        return std::make_unique<TriHexLattice>(cell_spacing, offset_x, offset_y);
    case ipMagmaHoneycomb:
        return std::make_unique<HexLattice>(cell_spacing, offset_x, offset_y, line_width);
    case ipMagmaTriangle:
    default:
        return std::make_unique<TriangleLattice>(cell_spacing, offset_x, offset_y);
    }
}

} // namespace magma
} // namespace Slic3r

#endif // slic3r_Magma_MagmaPatterns_hpp_
