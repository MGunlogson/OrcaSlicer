#include "ZoneInterior.hpp"

#include <libslic3r/OpenVDBUtils.hpp>
#include <libslic3r/TriangleMesh.hpp>
#include <libslic3r/SLA/Hollowing.hpp>
#include <libslic3r/Utils.hpp>

#include <openvdb/tools/LevelSetFilter.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/LevelSetUtil.h>   // sdfInteriorMask
#include <openvdb/tools/Morphology.h>      // dilateActiveValues
#include <openvdb/tools/Prune.h>           // pruneInactive

#include <boost/log/trivial.hpp>

#include <chrono>

namespace Slic3r {
namespace zone_boundary {

void smooth_interior(sla::Interior &interior, const TriangleMesh &original_mesh, int iterations)
{
    // Works on the grid only; interior.mesh may already be cleared by filter_thin_interior().
    if (!interior.gridptr || iterations <= 0)
        return;

    BOOST_LOG_TRIVIAL(info) << "Zone boundary: Applying constrained smoothing (" << iterations << " iterations)";

    // Clamping constraint: the original mesh, at the same voxel scale, offset inward by the shell thickness.
    float voxel_scale = float(interior.voxel_scale);
    float out_range = 3.0f;
    float in_range = float(interior.nb_in);

    auto original_grid = mesh_to_grid(original_mesh.its, {}, voxel_scale, out_range, in_range);
    if (!original_grid) {
        BOOST_LOG_TRIVIAL(warning) << "Zone boundary: Failed to create original mesh grid for smoothing constraint";
        return;
    }

    double thickness_offset = interior.thickness;
    auto valid_zone_grid = redistance_grid(*original_grid, -thickness_offset, in_range, in_range);
    if (!valid_zone_grid) {
        BOOST_LOG_TRIVIAL(warning) << "Zone boundary: Failed to create valid zone grid";
        return;
    }

    // Each iteration runs several smoothing passes, then clamps to the valid zone. Mean curvature
    // shrinks bumps inward and the clamp only blocks outward growth into the shell, so batching
    // passes loses nothing and saves CSG operations.
    const int smooth_passes_per_clamp = 5;

    // Convergence measure: sum of |SDF| over all active voxels. Cheap next to the smoothing passes.
    auto compute_l1_energy = [&]() -> double {
        double energy = 0;
        for (auto iter = interior.gridptr->cbeginValueOn(); iter; ++iter)
            energy += std::abs(iter.getValue());
        return energy;
    };

    double prev_energy = compute_l1_energy();
    double prev_rel_change = 0;

    BOOST_LOG_TRIVIAL(info) << "Zone boundary: Smoothing " << interior.gridptr->activeVoxelCount()
        << " active voxels, " << iterations << " max iterations x " << smooth_passes_per_clamp << " passes";

    auto t_smooth_start = std::chrono::high_resolution_clock::now();

    int actual_iterations = 0;
    for (int i = 0; i < iterations; ++i) {
        auto t_iter = std::chrono::high_resolution_clock::now();

        openvdb::tools::LevelSetFilter<openvdb::FloatGrid> filter(*interior.gridptr);

        for (int j = 0; j < smooth_passes_per_clamp; ++j) {
            filter.meanCurvature();
        }

        // csgIntersectionCopy leaves valid_zone_grid intact for the next iteration.
        interior.gridptr = openvdb::tools::csgIntersectionCopy(*interior.gridptr, *valid_zone_grid);

        ++actual_iterations;

        double curr_energy = compute_l1_energy();
        double rel_change = (prev_energy > 0) ? std::abs(curr_energy - prev_energy) / prev_energy : 0;

        auto iter_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - t_iter).count();

        BOOST_LOG_TRIVIAL(debug) << "Zone boundary: Smooth iter " << (i + 1) << "/" << iterations
            << " " << iter_ms << "ms, L1=" << curr_energy << ", rel_change=" << rel_change;

        // Stop once the relative change no longer halves per iteration: smoothing and clamping
        // have reached equilibrium. The first iteration is always a large change, so skip it.
        if (i > 0 && prev_rel_change > 0 && rel_change >= 0.5 * prev_rel_change) {
            BOOST_LOG_TRIVIAL(info) << "Zone boundary: Smooth plateaued at iteration " << (i + 1)
                << "/" << iterations << " (rel_change=" << rel_change << ")";
            break;
        }

        prev_rel_change = rel_change;
        prev_energy = curr_energy;
    }

    auto smooth_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - t_smooth_start).count();
    BOOST_LOG_TRIVIAL(info) << "Zone boundary: Smooth complete: " << actual_iterations
        << "/" << iterations << " iterations in " << smooth_ms << "ms";

    double adaptivity = 0.;
    interior.mesh = grid_to_mesh(*interior.gridptr, 0.0, adaptivity);

    if (!interior.mesh.empty())
        sla::postprocess_interior_mesh(interior.mesh);

    BOOST_LOG_TRIVIAL(info) << "Zone boundary: Constrained smoothing complete";
}

void filter_thin_interior(sla::Interior &interior, double min_width)
{
    if (!interior.gridptr || min_width <= 0)
        return;

    BOOST_LOG_TRIVIAL(info) << "Zone boundary: Filtering thin inner zone sections (min width: " << min_width << "mm)";

    // Inscribed sphere radius, in voxels
    float threshold = float(min_width / 2.0 * interior.voxel_scale);
    int threshold_voxels = int(std::ceil(threshold));

    // Thick core: voxels at least `threshold` inside the surface (SDF <= -threshold).
    auto thick_mask = openvdb::tools::sdfInteriorMask(*interior.gridptr, -threshold);

    if (!thick_mask || thick_mask->tree().activeVoxelCount() == 0) {
        BOOST_LOG_TRIVIAL(warning) << "Zone boundary: No regions thick enough to survive filtering (min_width=" << min_width << "mm)";
        interior.gridptr->clear();
        interior.mesh.clear();
        return;
    }

    // Dilate the core back out to the original surface.
    openvdb::tools::dilateActiveValues(thick_mask->tree(), threshold_voxels,
                                        openvdb::tools::NN_FACE_EDGE_VERTEX,
                                        openvdb::tools::PRESERVE_TILES);

    // Outside the mask, set voxels to background (outside); inside it keep the original SDF.
    float background = interior.gridptr->background();
    size_t removed_count = 0;
    size_t kept_count = 0;

    for (auto iter = interior.gridptr->beginValueOn(); iter; ++iter) {
        openvdb::Coord coord = iter.getCoord();
        if (!thick_mask->tree().isValueOn(coord)) {
            iter.setValue(background);
            removed_count++;
        } else {
            kept_count++;
        }
    }

    openvdb::tools::pruneInactive(interior.gridptr->tree());

    // Stale; smooth_interior() regenerates it.
    interior.mesh.clear();

    BOOST_LOG_TRIVIAL(info) << "Zone boundary: Thin inner zone filtering complete (kept " << kept_count << " voxels)";
}

bool debug_export_interior(const sla::Interior &interior, const std::string &stage_name, int object_id)
{
    if (interior.mesh.empty()) {
        BOOST_LOG_TRIVIAL(debug) << "Zone boundary debug: Skipping export of empty mesh for stage " << stage_name;
        return false;
    }

    std::string filename = debug_out_path("zone_shell_obj%d_%s.stl", object_id, stage_name.c_str());

    bool success = its_write_stl_ascii(filename.c_str(), "zone_interior", interior.mesh);

    if (success) {
        BOOST_LOG_TRIVIAL(info) << "Zone boundary debug: Exported interior mesh to " << filename
                                << " (" << interior.mesh.indices.size() << " triangles)";
    } else {
        BOOST_LOG_TRIVIAL(error) << "Zone boundary debug: Failed to export interior mesh to " << filename;
    }

    return success;
}

} // namespace zone_boundary
} // namespace Slic3r
