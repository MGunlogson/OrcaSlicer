#include "MagmaInterior.hpp"

#include <libslic3r/OpenVDBUtils.hpp>
#include <libslic3r/TriangleMesh.hpp>
#include <libslic3r/QuadricEdgeCollapse.hpp>
#include <libslic3r/Utils.hpp>

#include <openvdb/tools/LevelSetFilter.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/LevelSetUtil.h>   // sdfInteriorMask
#include <openvdb/tools/Morphology.h>      // dilateActiveValues
#include <openvdb/tools/Prune.h>           // pruneInactive

#include <boost/log/trivial.hpp>

namespace Slic3r {
namespace magma {

// Swap face winding to flip normals
static inline void swap_normals(indexed_triangle_set &its)
{
    for (auto &face : its.indices)
        std::swap(face(0), face(2));
}

void regenerate_mesh_from_grid(sla::Interior &interior)
{
    if (!interior.gridptr)
        return;

    double adaptivity = 0.;
    interior.mesh = grid_to_mesh(*interior.gridptr, 0.0, adaptivity);

    if (!interior.mesh.empty()) {
        swap_normals(interior.mesh);
    }

    interior.reset_accessor();
}

void smooth_interior(sla::Interior &interior, const TriangleMesh &original_mesh, int iterations)
{
    std::cerr << "Magma DEBUG smooth_interior: gridptr=" << (interior.gridptr ? "valid" : "null")
              << " iterations=" << iterations
              << " mesh_empty=" << interior.mesh.empty() << std::endl;

    // Note: interior.mesh may be empty if filter_thin_interior() was called first
    // We only need the grid for smoothing - mesh is regenerated at the end
    if (!interior.gridptr || iterations <= 0) {
        std::cerr << "Magma DEBUG smooth_interior: early return (no grid or iterations)" << std::endl;
        return;
    }

    BOOST_LOG_TRIVIAL(info) << "Magma: Applying constrained smoothing (" << iterations << " iterations)";

    // Create grid from original mesh at same voxel scale for clamping constraint.
    // The valid zone is the original mesh offset inward by shell thickness.
    float voxel_scale = float(interior.voxel_scale);
    float out_range = 3.0f;
    float in_range = float(interior.nb_in);

    auto original_grid = mesh_to_grid(original_mesh.its, {}, voxel_scale, out_range, in_range);
    if (!original_grid) {
        BOOST_LOG_TRIVIAL(warning) << "Magma: Failed to create original mesh grid for smoothing constraint";
        return;
    }

    // The valid zone grid: original mesh offset inward by thickness.
    // In SDF terms, we add the thickness offset to shift the zero level-set inward.
    // valid_zone represents the boundary that the interior must not cross.
    double thickness_offset = interior.thickness;
    auto valid_zone_grid = redistance_grid(*original_grid, -thickness_offset, in_range, in_range);
    if (!valid_zone_grid) {
        BOOST_LOG_TRIVIAL(warning) << "Magma: Failed to create valid zone grid";
        return;
    }

    // Apply constrained mean curvature smoothing iterations.
    // Each iteration: multiple smooth passes, then clamp to valid zone.
    // Batching smooth passes before clamping is more effective because:
    // - Mean curvature shrinks bumps INWARD (convex -> shrink)
    // - Clamping only prevents OUTWARD expansion (into shell zone)
    // - Clamping does NOT block bump removal
    // - Batching = more effective smoothing per CSG operation
    //
    // We use csgIntersectionCopy which:
    // - Takes const references (doesn't destroy valid_zone_grid)
    // - Maintains proper level set semantics (signedFloodFill, narrow band)
    // - Is internally parallelized with TBB

    const int smooth_passes_per_clamp = 5;  // Batch 5 smooth passes before clamping

    for (int i = 0; i < iterations; ++i) {
        openvdb::tools::LevelSetFilter<openvdb::FloatGrid> filter(*interior.gridptr);

        // Apply multiple smoothing passes before clamping
        for (int j = 0; j < smooth_passes_per_clamp; ++j) {
            filter.meanCurvature();
        }

        // Clamp to valid zone using csgIntersectionCopy (preserves both inputs)
        interior.gridptr = openvdb::tools::csgIntersectionCopy(*interior.gridptr, *valid_zone_grid);
    }

    std::cerr << "Magma DEBUG smooth_interior: applied " << (iterations * smooth_passes_per_clamp)
              << " smoothing passes (" << iterations << " iterations x " << smooth_passes_per_clamp << " passes)" << std::endl;

    // Convert smoothed grid back to mesh
    double adaptivity = 0.;
    interior.mesh = grid_to_mesh(*interior.gridptr, 0.0, adaptivity);

    std::cerr << "Magma DEBUG smooth_interior: after grid_to_mesh, mesh_empty=" << interior.mesh.empty()
              << " triangles=" << interior.mesh.indices.size() << std::endl;

    if (!interior.mesh.empty()) {
        // Post-process mesh as in generate_interior
        swap_normals(interior.mesh);
        float loss_less_max_error = 2 * std::numeric_limits<float>::epsilon();
        its_quadric_edge_collapse(interior.mesh, 0U, &loss_less_max_error);
        its_compactify_vertices(interior.mesh);
        its_merge_vertices(interior.mesh);
        swap_normals(interior.mesh);
    }

    std::cerr << "Magma DEBUG smooth_interior: final mesh_empty=" << interior.mesh.empty()
              << " triangles=" << interior.mesh.indices.size() << std::endl;

    BOOST_LOG_TRIVIAL(info) << "Magma: Constrained smoothing complete";
}

void filter_thin_interior(sla::Interior &interior, double min_width)
{
    if (!interior.gridptr || min_width <= 0)
        return;

    BOOST_LOG_TRIVIAL(info) << "Magma: Filtering thin yolk sections (min width: " << min_width << "mm)";

    // Convert min_width to voxel scale
    // The threshold is half the min_width (radius of inscribed sphere)
    float threshold = float(min_width / 2.0 * interior.voxel_scale);
    int threshold_voxels = int(std::ceil(threshold));

    std::cerr << "Magma DEBUG filter_thin_interior: min_width=" << min_width
              << "mm threshold=" << threshold << " voxels (" << threshold_voxels << " int)"
              << " voxel_scale=" << interior.voxel_scale << std::endl;

    // 1. Extract mask of "thick core" - voxels where inscribed sphere radius >= threshold
    //    sdfInteriorMask returns mask where SDF <= isovalue
    //    Negative threshold = voxels at least threshold from boundary
    auto thick_mask = openvdb::tools::sdfInteriorMask(*interior.gridptr, -threshold);

    if (!thick_mask || thick_mask->tree().activeVoxelCount() == 0) {
        // No thick regions survive - clear the interior
        BOOST_LOG_TRIVIAL(warning) << "Magma: No regions thick enough to survive filtering (min_width=" << min_width << "mm)";
        std::cerr << "Magma DEBUG filter_thin_interior: no thick regions survive, clearing interior" << std::endl;
        interior.gridptr->clear();
        interior.mesh.clear();
        return;
    }

    size_t core_voxels = thick_mask->tree().activeVoxelCount();
    std::cerr << "Magma DEBUG filter_thin_interior: thick core has " << core_voxels << " voxels" << std::endl;

    // 2. Dilate mask back to reach original surface (BINARY dilation = EXACT)
    openvdb::tools::dilateActiveValues(thick_mask->tree(), threshold_voxels,
                                        openvdb::tools::NN_FACE_EDGE_VERTEX,
                                        openvdb::tools::PRESERVE_TILES);

    size_t dilated_voxels = thick_mask->tree().activeVoxelCount();
    std::cerr << "Magma DEBUG filter_thin_interior: dilated mask has " << dilated_voxels << " voxels" << std::endl;

    // 3. Use mask to carve original level set
    //    Keep original SDF where mask is active, set to background (outside) elsewhere
    float background = interior.gridptr->background();
    size_t removed_count = 0;
    size_t kept_count = 0;

    for (auto iter = interior.gridptr->beginValueOn(); iter; ++iter) {
        openvdb::Coord coord = iter.getCoord();
        if (!thick_mask->tree().isValueOn(coord)) {
            iter.setValue(background);  // Set to outside
            removed_count++;
        } else {
            kept_count++;
        }
    }

    std::cerr << "Magma DEBUG filter_thin_interior: kept " << kept_count
              << " voxels, removed " << removed_count << " voxels" << std::endl;

    // Prune inactive voxels
    openvdb::tools::pruneInactive(interior.gridptr->tree());

    // Clear stale mesh - smooth_interior() will regenerate it
    interior.mesh.clear();

    BOOST_LOG_TRIVIAL(info) << "Magma: Thin yolk filtering complete (kept " << kept_count << " voxels)";
}

bool debug_export_interior(const sla::Interior &interior, const std::string &stage_name, int object_id)
{
    if (interior.mesh.empty()) {
        BOOST_LOG_TRIVIAL(debug) << "Magma debug: Skipping export of empty mesh for stage " << stage_name;
        return false;
    }

    std::string filename = debug_out_path("magma_shell_obj%d_%s.stl", object_id, stage_name.c_str());

    bool success = its_write_stl_ascii(filename.c_str(), "magma_interior", interior.mesh);

    if (success) {
        BOOST_LOG_TRIVIAL(info) << "Magma debug: Exported interior mesh to " << filename
                                << " (" << interior.mesh.indices.size() << " triangles)";
    } else {
        BOOST_LOG_TRIVIAL(error) << "Magma debug: Failed to export interior mesh to " << filename;
    }

    return success;
}

} // namespace magma
} // namespace Slic3r
