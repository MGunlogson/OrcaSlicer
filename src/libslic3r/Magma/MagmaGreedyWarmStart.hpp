#ifndef slic3r_Magma_MagmaGreedyWarmStart_hpp_
#define slic3r_Magma_MagmaGreedyWarmStart_hpp_

#include "MagmaTriangleCell.hpp"
#include "MagmaTubeMap.hpp"
#include "MagmaTubeSolver.hpp"

#include <cstdint>
#include <unordered_map>
#include <vector>

namespace Slic3r {
namespace magma {

/// Greedy warm start: assign tubes using a most-constrained-first heuristic.
/// Populates committed (must be pre-sized to edges.size()) with layer-aligned
/// CommittedSegments that respect min/max height, run boundaries, and per-cell
/// NoOverlap. Fast (100-300ms for 1000+ cells) and deterministic.
void greedy_warm_start(
    const std::unordered_map<TriangleCell, CellPresence, TriangleCellHash> &cells,
    const std::vector<EdgeData>                                             &edges,
    const std::unordered_map<TriangleCell, std::vector<size_t>, TriangleCellHash> &cell_edges,
    const MicronTables                                                      &um,
    int64_t min_h_um,
    int64_t max_h_um,
    std::vector<std::vector<CommittedSegment>>                              &committed);

} // namespace magma
} // namespace Slic3r

#endif // slic3r_Magma_MagmaGreedyWarmStart_hpp_
