#ifndef slic3r_Magma_MagmaInterior_hpp_
#define slic3r_Magma_MagmaInterior_hpp_

#include <string>
#include <libslic3r/SLA/Interior.hpp>

namespace Slic3r {

class TriangleMesh;

namespace magma {

// Regenerate mesh from grid. Call this after filter_thin_interior() if you need
// to capture the mesh state before smoothing.
void regenerate_mesh_from_grid(sla::Interior &interior);

// Apply constrained mean curvature smoothing to the interior boundary.
// Smooths stair-step artifacts while ensuring shell thickness is maintained.
// original_mesh: the original mesh (for computing thickness constraint)
// iterations: number of smoothing iterations (default 5)
void smooth_interior(sla::Interior &interior, const TriangleMesh &original_mesh, int iterations = 5);

// Filter out thin yolk sections using morphological reconstruction.
// Removes regions where the interior is thinner than min_width in any direction.
// This eliminates small disconnected islands and thin protrusions that would
// create unusable infill zones.
// Algorithm: erode to find thick core, over-dilate, intersect with original
// to preserve exact surface detail in thick regions.
// min_width: minimum thickness in mm (0 to disable)
void filter_thin_interior(sla::Interior &interior, double min_width);

// Export interior mesh to STL for debugging.
// stage_name: descriptive name for this processing stage (e.g., "1_initial", "2_filtered")
// object_id: ID of the PrintObject for unique filenames
// Returns true if export was successful, false if mesh is empty or export failed.
bool debug_export_interior(const sla::Interior &interior, const std::string &stage_name, int object_id = 0);

} // namespace magma
} // namespace Slic3r

#endif // slic3r_Magma_MagmaInterior_hpp_
