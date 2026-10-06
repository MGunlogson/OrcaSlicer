#ifndef slic3r_ZoneBoundary_ZoneInterior_hpp_
#define slic3r_ZoneBoundary_ZoneInterior_hpp_

#include <string>
#include <libslic3r/SLA/Interior.hpp>

namespace Slic3r {

class TriangleMesh;

namespace zone_boundary {

// Mean curvature smoothing of the interior boundary, clamped so the shell never gets thinner than
// interior.thickness from original_mesh. Runs at most `iterations` rounds, stopping once converged.
// Regenerates interior.mesh from the grid.
void smooth_interior(sla::Interior &interior, const TriangleMesh &original_mesh, int iterations = 5);

// Remove interior regions thinner than min_width (mm, 0 = off): islands and protrusions too thin
// to be a usable inner zone. Erodes to the thick core, dilates it back and keeps the original SDF
// inside that mask, so thick regions keep their exact surface. Works on the grid only and clears
// interior.mesh.
void filter_thin_interior(sla::Interior &interior, double min_width);

// Debug: export interior.mesh to STL. False if the mesh is empty or the write failed.
bool debug_export_interior(const sla::Interior &interior, const std::string &stage_name, int object_id = 0);

} // namespace zone_boundary
} // namespace Slic3r

#endif // slic3r_ZoneBoundary_ZoneInterior_hpp_
