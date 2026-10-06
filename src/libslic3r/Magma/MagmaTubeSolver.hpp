#ifndef slic3r_Magma_MagmaTubeSolver_hpp_
#define slic3r_Magma_MagmaTubeSolver_hpp_

#include "MagmaTriangleCell.hpp"
#include "MagmaTubeMap.hpp"

#include <cstdint>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace Slic3r {
namespace magma {

// ============================================================================
// CellEdge — canonical adjacent cell pair (a < b)
// ============================================================================

struct CellEdge {
    CellId a, b;

    CellEdge() = default;
    CellEdge(const CellId &x, const CellId &y)
        : a(std::min(x, y)), b(std::max(x, y)) {}

    bool operator==(const CellEdge &o) const { return a == o.a && b == o.b; }
    bool operator<(const CellEdge &o) const {
        if (a != o.a) return a < o.a;
        return b < o.b;
    }
};

struct CellEdgeHash {
    size_t operator()(const CellEdge &e) const {
        CellIdHash h;
        size_t seed = h(e.a);
        seed ^= h(e.b) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

// ============================================================================
// Run — contiguous shared-presence range in microns
// ============================================================================

struct Run {
    int     start_layer, end_layer; // inclusive layer range
    int64_t start_um, end_um;       // bottom_um[start_layer], top_um[end_layer]
};

// ============================================================================
// EdgeData — pre-computed per-edge data
// ============================================================================

struct EdgeData {
    CellEdge         edge;
    std::vector<Run> runs; // contiguous shared-presence ranges (split where either cell's presence breaks)
};

// ============================================================================
// CommittedSegment — a tube assignment in micron space
// ============================================================================

struct CommittedSegment {
    int64_t start_um, end_um;
};

// ============================================================================
// MicronTables — layer boundary Z values in integer microns
// ============================================================================
// Contiguous by construction: bottom_um[L+1] == top_um[L].

struct MicronTables {
    std::vector<int64_t> top_um;    // top_um[L] = llround(print_z * 1000)
    std::vector<int64_t> bottom_um; // bottom_um[0] from bottom_z; L>0: top_um[L-1]

    std::unordered_map<int64_t, int> bottom_to_layer; // reverse lookup
    std::unordered_map<int64_t, int> top_to_layer;    // reverse lookup
};

// ============================================================================
// Block — 3D region of cells x layers, solved as one CP-SAT model
// ============================================================================

struct Block {
    std::vector<size_t>                                    edge_indices; // into m_edges
    std::unordered_set<CellId, CellIdHash>     cells;
    int     z_start_layer, z_end_layer;
    int64_t z_start_um, z_end_um;
};

// ============================================================================
// BlockResult — output of solving one block
// ============================================================================

// NoWork and Failed both leave the committed state alone; only Failed is reported to the user.
enum class BlockOutcome {
    Solved,  // commit these segments in place of what is there
    NoWork,  // nothing to decide (no tubes in range, or no pairs to separate) -- keep as-is
    Failed,  // solver returned no solution (INFEASIBLE, UNKNOWN, ...) -- keep as-is, and say so
};

struct BlockResult {
    BlockOutcome outcome = BlockOutcome::Failed;
    std::vector<std::pair<size_t, CommittedSegment>> segments; // (edge_idx, segment)
};

// ============================================================================
// ValidationResult — output of MagmaTubeSolver::validate_and_report()
// ============================================================================

struct ValidationResult {
    int bad_short = 0, bad_long = 0, bad_range = 0, bad_presence = 0;
    int overlap_cell = 0, overlap_edge = 0;
    int total_issues() const {
        return bad_range + bad_short + bad_long +
               bad_presence + overlap_cell + overlap_edge;
    }
};

// Smallest positive layer height at or above `first_layer`, or 0.0 when there is none.
// layer_data is indexed by absolute Layer::id(), so with a raft the rows below first_layer
// are zero-height placeholders and must not seed the minimum.
double min_positive_layer_height(const std::vector<LayerData> &layer_data, int first_layer);

// ============================================================================
// MagmaTubeSolver — CP-SAT interval scheduling solver for tube assignment
// ============================================================================

class MagmaTubeSolver {
public:
    using ProgressFn      = std::function<void(int, int)>; // (current, total)
    using ThrowIfCanceled = std::function<void()>;

    MagmaTubeSolver(
        const MagmaLattice &lattice,
        const std::unordered_map<CellId, CellPresence, CellIdHash> &cells,
        const std::vector<LayerData> &layer_data,
        int    first_layer,
        double min_tube_height_mm,
        double max_tube_height_mm,
        int    num_layers,
        MagmaTubeSolverMode mode = MagmaTubeSolverMode::CoverageStagger,
        double solver_timeout_sec = 20.0);

    /// Run the solver. Populates out_pairs and out_cell_pair_index.
    void solve(
        std::vector<UTubePair> &out_pairs,
        std::unordered_map<CellId, std::vector<int>, CellIdHash> &out_cell_pair_index,
        ProgressFn progress_fn = nullptr,
        ThrowIfCanceled throw_if_canceled = nullptr);

    /// Blocks whose solve found no solution; their committed state is kept.
    int failed_block_count() const { return m_failed_blocks; }

private:
    // Pre-computation
    void build_micron_tables();
    void build_edges();
    void build_blocks(int off_z, std::vector<Block> &out) const;

    // One weighted objective per block: coverage dominates; tube count and stagger break ties.
    BlockResult solve_block(const Block &block) const;
    void        solve_pass(const std::vector<Block> &blocks,
                           std::function<void(int)> block_done_fn = nullptr);
    void        commit_result(const Block &block, const BlockResult &result);

    // Output conversion
    void extract_results(
        std::vector<UTubePair> &out_pairs,
        std::unordered_map<CellId, std::vector<int>,
                           CellIdHash> &out_index) const;

    // Diagnostics, run under label GREEDY after the warm start and CPSAT after the solve, so
    // the two compare line for line.

    /// Validate committed segments (height bounds, presence, per-cell and per-edge overlap),
    /// then report coverage, tube-length distribution and stagger quality.
    void validate_and_report(const char *label) const;

    /// Fraction of each cell's presence range that ended up inside a tube.
    void report_coverage(
        const std::unordered_map<CellId, std::vector<std::pair<int, int>>,
                                 CellIdHash> &cell_segments,
        const char *label) const;

    /// Tube-height distribution and same-edge abutment count.
    void report_lengths(const char *label) const;

    /// Stagger quality: how far apart tube boundaries sit within each cell's
    /// Ring-1 neighbourhood, and how many cells share the worst Z plane.
    void report_stagger(const char *label) const;

    // Input (caller owns the data).
    // Lattice topology is offset-independent, so any layer's lattice serves.
    const MagmaLattice &m_lattice;
    const std::unordered_map<CellId, CellPresence, CellIdHash> &m_cells;
    const std::vector<LayerData> &m_layer_data;

    // Config
    double m_min_h_mm;
    double m_max_h_mm;
    int    m_num_layers;
    // First m_layer_data row backed by a real layer (non-zero with a raft, since Layer::id()
    // is absolute). Rows below it are placeholders and must not be read.
    int    m_first_layer;
    int    m_z_window; // Z block size in layers
    MagmaTubeSolverMode m_mode;
    double m_timeout_sec;

    // Pre-computed
    MicronTables                                            m_um;
    std::vector<EdgeData>                                   m_edges;
    // Cell → edge indices involving that cell (reverse lookup for stagger)
    std::unordered_map<CellId, std::vector<size_t>, CellIdHash> m_cell_edges;

    // Committed assignments (updated between passes)
    std::vector<std::vector<CommittedSegment>> m_committed; // indexed by edge idx
    int m_failed_blocks = 0; // BlockOutcome::Failed count

    // Computed at solve time.
    // The time budget is re-divided before every block (budget left / blocks left), so time
    // unused by easy blocks passes to the hard ones.
    double m_per_block_timeout = 10.0;  // this block's slice, set by solve_pass
    double m_budget_left_sec   = 0.0;   // unspent budget
    int    m_blocks_left       = 0;     // blocks still to solve, across all Z levels
    int    m_cpsat_workers = 8;         // set to TBB max concurrency by solve()

    // Constants
    static constexpr int    R              = 16;  // XY block size
    static constexpr int    R_OVERLAP      = 2;   // XY overlap between adjacent blocks
};

} // namespace magma
} // namespace Slic3r

#endif // slic3r_Magma_MagmaTubeSolver_hpp_
