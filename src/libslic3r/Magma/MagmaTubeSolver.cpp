#include "MagmaTubeSolver.hpp"
#include "MagmaGreedyWarmStart.hpp"

#include <boost/log/trivial.hpp>
#include <tbb/task_arena.h> // for max_concurrency()

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <array>
#include <map>
#include <set>

// OR-Tools CP-SAT solver
#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/sat/cp_model_checker.h"

namespace Slic3r {
namespace magma {

// Alias rather than a using-directive: operations_research::sat::Model would clash with
// Slic3r::Model.
namespace ortsat = operations_research::sat;

double min_positive_layer_height(const std::vector<LayerData> &layer_data, int first_layer)
{
    double min_lh = 0.0;
    for (int L = std::max(0, first_layer); L < int(layer_data.size()); ++L)
        if (layer_data[L].height > 0 && (min_lh <= 0.0 || layer_data[L].height < min_lh))
            min_lh = layer_data[L].height;
    return min_lh;
}

// ============================================================================
// Constructor
// ============================================================================

MagmaTubeSolver::MagmaTubeSolver(
    const MagmaLattice &lattice,
    const std::unordered_map<CellId, CellPresence, CellIdHash> &cells,
    const std::vector<LayerData> &layer_data,
    int    first_layer,
    double min_tube_height_mm,
    double max_tube_height_mm,
    int    num_layers,
    MagmaTubeSolverMode mode,
    double solver_timeout_sec)
    : m_lattice(lattice)
    , m_cells(cells)
    , m_layer_data(layer_data)
    , m_min_h_mm(min_tube_height_mm)
    , m_max_h_mm(max_tube_height_mm)
    , m_num_layers(num_layers)
    , m_first_layer(std::max(0, first_layer))
    , m_z_window(0)
    , m_mode(mode)
    , m_timeout_sec(solver_timeout_sec)
{}

// ============================================================================
// solve — greedy warm start, then one CP-SAT sweep over all blocks
// ============================================================================

void MagmaTubeSolver::solve(
    std::vector<UTubePair> &out_pairs,
    std::unordered_map<CellId, std::vector<int>, CellIdHash> &out_cell_pair_index,
    ProgressFn progress_fn,
    ThrowIfCanceled throw_if_canceled)
{
    auto t_start = std::chrono::high_resolution_clock::now();

    // Smallest layer height, for a conservative layer count per tube. Only non-positive with no
    // layers at all, and then there is nothing to assign.
    const double min_lh = min_positive_layer_height(m_layer_data, m_first_layer);
    if (min_lh <= 0.0) {
        BOOST_LOG_TRIVIAL(error) << "Magma solver: no layer with a positive height; "
                                    "skipping tube assignment.";
        return;
    }

    build_micron_tables();
    build_edges();
    if (throw_if_canceled) throw_if_canceled();

    // Greedy fills m_committed; CP-SAT then refines it, using it as a solution hint.
    {
        const int64_t min_h_um = llround(m_min_h_mm * 1000.0);
        const int64_t max_h_um = llround(m_max_h_mm * 1000.0);
        greedy_warm_start(m_lattice, m_cells, m_edges, m_cell_edges, m_um,
                          min_h_um, max_h_um, m_committed);
        validate_and_report("GREEDY");
    }
    if (throw_if_canceled) throw_if_canceled();

    // Z window and stride, derived from max tube height:
    //   window  = 4 × max_h_layers  — room for 3-4 stacked tubes
    //   overlap = 2 × max_h_layers  — every tube fully in ≥2 Z levels
    //   stride  = 2 × max_h_layers  — window minus overlap
    int z_stride;
    {
        int max_h_layers = std::max(1, static_cast<int>(std::ceil(m_max_h_mm / min_lh)));
        m_z_window = 4 * max_h_layers;
        z_stride   = std::max(1, 2 * max_h_layers);
    }

    // Z extent of actual edge data (not all layers may have cells)
    int max_edge_layer = 0;
    for (const auto &[cell, presence] : m_cells)
        max_edge_layer = std::max(max_edge_layer, presence.last_layer);

    int num_z_levels = 0;
    for (int z = 0; z <= max_edge_layer; z += z_stride)
        ++num_z_levels;

    // Give CP-SAT all available cores (blocks run sequentially).
    m_cpsat_workers = tbb::this_task_arena::max_concurrency();

    // Blocks depend only on immutable state, so build them once up front; the total also
    // sizes the time budget.
    std::vector<std::vector<Block>> z_levels;
    int total_blocks = 0;
    for (int z_off = 0; z_off <= max_edge_layer; z_off += z_stride) {
        std::vector<Block> level;
        build_blocks(z_off, level);
        total_blocks += static_cast<int>(level.size());
        z_levels.push_back(std::move(level));
    }
    m_budget_left_sec   = m_timeout_sec;
    m_blocks_left       = total_blocks;
    m_per_block_timeout = m_timeout_sec / std::max(1, total_blocks);

    BOOST_LOG_TRIVIAL(info) << "MagmaTubeSolver: " << m_edges.size() << " edges"
        << ", Z_window=" << m_z_window << " layers"
        << ", Z_stride=" << z_stride
        << ", R=" << R
        << ", " << num_z_levels << " Z levels"
        << ", " << total_blocks << " total blocks"
        << ", budget=" << std::fixed << std::setprecision(1) << m_timeout_sec << "s"
        << " (redivided per block; " << m_per_block_timeout << "s if spread evenly)"
        << ", cpsat_workers=" << m_cpsat_workers
        << ", max_layer=" << max_edge_layer;

    // One XY pass per Z level over overlapping blocks, solved sequentially with CP-SAT using
    // all cores within each. Both modes run it; the mode only decides whether the stagger
    // term is in the objective.
    {
        int work_done = 0;
        auto step = [&](int) {
            ++work_done;
            if (progress_fn) progress_fn(work_done, total_blocks);
            if (throw_if_canceled) throw_if_canceled();
        };

        if (progress_fn) progress_fn(0, total_blocks);
        for (const std::vector<Block> &level : z_levels)
            solve_pass(level, step);

        validate_and_report("CPSAT");
    }

    extract_results(out_pairs, out_cell_pair_index);

    auto t_end = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();
    BOOST_LOG_TRIVIAL(info) << "MagmaTubeSolver: "
        << (m_mode == MagmaTubeSolverMode::CoverageStagger ? "coverage+stagger" : "coverage")
        << " in " << ms << "ms, " << out_pairs.size() << " pairs placed";
}

// ============================================================================
// build_micron_tables
// ============================================================================

void MagmaTubeSolver::build_micron_tables()
{
    int n = m_num_layers;
    m_um.top_um.resize(n);
    m_um.bottom_um.resize(n);

    // Only rows from m_first_layer up are real; placeholder rows would all map Z=0 in the
    // reverse maps below.
    for (int L = m_first_layer; L < n; ++L)
        m_um.top_um[L] = llround(m_layer_data[L].print_z * 1000.0);

    // The first real layer's bottom is its own bottom_z (the raft top, with a raft).
    if (m_first_layer < n)
        m_um.bottom_um[m_first_layer] = llround(m_layer_data[m_first_layer].bottom_z() * 1000.0);
    for (int L = m_first_layer + 1; L < n; ++L)
        m_um.bottom_um[L] = m_um.top_um[L - 1];

    // Reverse maps (same integers — guaranteed to match)
    for (int L = m_first_layer; L < n; ++L) {
        m_um.bottom_to_layer[m_um.bottom_um[L]] = L;
        m_um.top_to_layer[m_um.top_um[L]] = L;
    }
}

// ============================================================================
// build_edges
// ============================================================================

void MagmaTubeSolver::build_edges()
{
    const int64_t min_h_um = llround(m_min_h_mm * 1000.0);

    std::unordered_set<CellEdge, CellEdgeHash> seen;

    for (const auto &[cell, presence] : m_cells) {
        for (const CellId &nbr : m_lattice.neighbors(cell)) {
            if (m_cells.find(nbr) == m_cells.end())
                continue;
            CellEdge edge(cell, nbr);
            if (!seen.insert(edge).second)
                continue;

            EdgeData ed;
            ed.edge = edge;

            const CellPresence &pa = m_cells.at(edge.a);
            const CellPresence &pb = m_cells.at(edge.b);
            int range_start = std::max(pa.first_layer, pb.first_layer);
            int range_end   = std::min(pa.last_layer, pb.last_layer);

            // Walk shared presence, split into contiguous runs
            Run current{-1, -1, 0, 0};
            auto flush_run = [&]() {
                if (current.start_layer < 0)
                    return;
                current.start_um = m_um.bottom_um[current.start_layer];
                current.end_um   = m_um.top_um[current.end_layer];
                if (current.end_um - current.start_um >= min_h_um)
                    ed.runs.push_back(current);
                current = {-1, -1, 0, 0};
            };

            for (int i = range_start; i <= range_end; ++i) {
                if (pa.present(i) && pb.present(i)) {
                    if (current.start_layer < 0)
                        current.start_layer = i;
                    current.end_layer = i;
                } else {
                    flush_run();
                }
            }
            flush_run();

            if (!ed.runs.empty())
                m_edges.push_back(std::move(ed));
        }
    }

    // Initialize committed segments (empty for all edges)
    m_committed.resize(m_edges.size());

    // Build reverse lookup: cell → edge indices
    for (size_t ei = 0; ei < m_edges.size(); ++ei) {
        m_cell_edges[m_edges[ei].edge.a].push_back(ei);
        m_cell_edges[m_edges[ei].edge.b].push_back(ei);
    }

    BOOST_LOG_TRIVIAL(debug) << "MagmaTubeSolver: built " << m_edges.size() << " edges";
}

// ============================================================================
// build_blocks
// ============================================================================

void MagmaTubeSolver::build_blocks(int off_z, std::vector<Block> &out) const
{
    out.clear();

    // Assign each cell to an XY block
    struct XYKey {
        int bx, by;
        bool operator==(const XYKey &o) const { return bx == o.bx && by == o.by; }
    };
    struct XYKeyHash {
        size_t operator()(const XYKey &k) const {
            return std::hash<int>()(k.bx) ^ (std::hash<int>()(k.by) << 16);
        }
    };

    // Assign each cell to its primary block, plus adjacent blocks if it falls
    // in the R_OVERLAP overlap zone. Stride = R - R_OVERLAP.
    const int stride = std::max(1, R - R_OVERLAP);

    auto floor_div = [](int a, int s) { return (a >= 0) ? a / s : (a - s + 1) / s; };

    // Anchor the grid on the part's minimum cell, not the lattice origin, so blocking is
    // invariant under moving the object and a part up to R cells across is one block.
    int min_a = std::numeric_limits<int>::max();
    int min_b = std::numeric_limits<int>::max();
    for (const auto &[cell, _] : m_cells) {
        min_a = std::min(min_a, cell.a);
        min_b = std::min(min_b, cell.b);
    }

    std::unordered_map<XYKey, std::vector<CellId>, XYKeyHash> xy_groups;
    for (const auto &[cell, _] : m_cells) {
        const int sa = cell.a - min_a;
        const int sb = cell.b - min_b;
        xy_groups[{floor_div(sa, stride), floor_div(sb, stride)}].push_back(cell);
    }

    // Overlap, in a second pass so it only widens blocks that already hold cells of their own.
    // A block made only of a neighbour's margin would be a subset of that neighbour, and
    // re-solving it could only repeat or undo the neighbour's work.
    for (const auto &[cell, _] : m_cells) {
        const int sa = cell.a - min_a, sb = cell.b - min_b;
        const int bx = floor_div(sa, stride), by = floor_div(sb, stride);
        // Block bx-1 spans [(bx-1)*stride, (bx-1)*stride + R) = up to bx*stride + R_OVERLAP,
        // so a cell also belongs to bx-1 exactly when its residual is below R_OVERLAP.
        const bool overlap_a = (sa - bx * stride) < R_OVERLAP;
        const bool overlap_b = (sb - by * stride) < R_OVERLAP;
        auto extend = [&](XYKey k) {
            auto it = xy_groups.find(k);
            if (it != xy_groups.end()) it->second.push_back(cell);
        };
        if (overlap_a)              extend({bx - 1, by});
        if (overlap_b)              extend({bx, by - 1});
        if (overlap_a && overlap_b) extend({bx - 1, by - 1});
    }

    // Single Z slice: off_z to off_z + z_window - 1
    // (Z levels are iterated by the caller in solve())
    int z_start = std::max(m_first_layer, off_z);   // never descend into placeholder rows
    int z_end   = std::min(off_z + m_z_window - 1, m_num_layers - 1);
    if (z_start >= m_num_layers || z_start > z_end) return;

    int64_t z_start_um = m_um.bottom_um[z_start];
    int64_t z_end_um   = m_um.top_um[z_end];

    // Deterministic order, low corner first: each block is solved against what its
    // predecessors committed, so order affects the result.
    std::vector<XYKey> keys;
    keys.reserve(xy_groups.size());
    for (const auto &[k, _] : xy_groups) keys.push_back(k);
    std::sort(keys.begin(), keys.end(), [](const XYKey &l, const XYKey &r) {
        return l.by != r.by ? l.by < r.by : l.bx < r.bx;
    });

    // Build blocks: one per XY group at this Z slice
    for (const XYKey &xy_key : keys) {
        const std::vector<CellId> &cells = xy_groups.at(xy_key);
        std::unordered_set<CellId, CellIdHash> cell_set(cells.begin(), cells.end());

        Block block;
        block.cells = cell_set;
        block.z_start_layer = z_start;
        block.z_end_layer   = z_end;
        block.z_start_um    = z_start_um;
        block.z_end_um      = z_end_um;

        // Collect edges via cell reverse lookup (avoids scanning all edges)
        std::unordered_set<size_t> seen_edges;
        for (const CellId &cell : cells) {
            auto it = m_cell_edges.find(cell);
            if (it == m_cell_edges.end()) continue;
            for (size_t ei : it->second) {
                if (!seen_edges.insert(ei).second) continue;
                const EdgeData &ed = m_edges[ei];
                if (cell_set.count(ed.edge.a) == 0 || cell_set.count(ed.edge.b) == 0)
                    continue;

                bool overlaps_z = false;
                for (const Run &r : ed.runs) {
                    if (r.end_um > z_start_um && r.start_um < z_end_um) {
                        overlaps_z = true;
                        break;
                    }
                }
                if (overlaps_z)
                    block.edge_indices.push_back(ei);
            }
        }

        if (block.edge_indices.empty()) continue;
        // Edge order is variable order in the model, which affects CP-SAT's search.
        std::sort(block.edge_indices.begin(), block.edge_indices.end());
        out.push_back(std::move(block));
    }

    BOOST_LOG_TRIVIAL(debug) << "MagmaTubeSolver: Z level " << off_z
        << " -> " << out.size() << " blocks";
}

// ============================================================================
// SegVars — one candidate tube slot in the CP-SAT model
// ============================================================================

struct SegVars {
    ortsat::BoolVar     active;
    ortsat::IntVar      start;
    ortsat::IntVar      end;
    ortsat::IntVar      size;
    size_t      edge_idx;
    int64_t     run_start_um;  // which run this slot belongs to
    int64_t     run_end_um;
    int         slot_idx;      // k within the run (0-based)
    operations_research::Domain contrib_dom; // {0} u feasible sizes
};

namespace {

// Distance in microns beyond which two tube boundaries no longer count as close. A modelling
// bound, not a setting: it only sets where the "closer is worse" penalty reaches zero.
//
// Must be at most half the minimum tube height. Two neighbouring columns with boundaries
// spaced s apart and offset by d give each boundary partners at d and s-d; if a horizon h
// reaches both, their penalties sum to 2h-s regardless of d, and the objective goes flat.
// Since s >= min_h, h = min_h/2 keeps only the nearer partner in range.
inline int64_t stagger_horizon_um(int64_t min_h_um) { return min_h_um / 2; }

// Where slot k's boundaries can land. Activation is a prefix (see the symmetry breaking in
// solve_block), so an active slot k has k tubes of at least min_h below it in the run.
// There is no tighter upper bound than the run end: slots above k need not be active, and
// excluding a legal position here makes the model infeasible.
struct SlotRange { int64_t lo, hi; };

SlotRange boundary_range(const SegVars &sv, bool is_end, int64_t min_h_um)
{
    const int64_t below = int64_t(sv.slot_idx) + (is_end ? 1 : 0);
    int64_t lo = sv.run_start_um + below * min_h_um;
    const int64_t hi = sv.run_end_um;
    if (lo > hi) lo = sv.run_start_um;
    return {lo, hi};
}

// How close to optimal the stagger term must be proven before the solver may stop, in
// fully-coincident boundary pairs. W_COVERAGE is inflated by the same amount, so raising it
// trades stagger quality for solve time without ever costing coverage.
constexpr int64_t STAGGER_GAP_PAIRS = 1;

// Crowding levels scored per window. J levels reproduce C(N,2) exactly up to N = J+1 and grow
// linearly above it, which is still monotone.
constexpr int STAGGER_HINGE_LEVELS = 4;

// Stagger penalty: slide a horizon-tall window up the layer grid and charge crowding.
//
// Equivalent to the per-pair tent max(0, horizon - |b1 - b2|): measured in grid steps, two
// boundaries d apart share (horizon - d) windows, so
//
//     sum over windows of C(N, 2)  ==  sum over pairs of max(0, horizon - d)
//
// and C(N,2) = sum over j>=1 of max(0, N - j), each hinge one variable and one linear row.
// Counting with convex hinges relaxes far better than the concave |b1 - b2|, on which the
// solver has to branch for every pair.
//
// The windows overlap, one per grid position; edge-to-edge tiling would let two boundaries a
// layer apart straddle a seam and score nothing. Frozen boundaries just add a constant to a
// position's count.
ortsat::LinearExpr build_stagger_penalty(
    ortsat::CpModelBuilder                                                       &model,
    const Block                                                          &block,
    const std::vector<SegVars>                                           &segs,
    const std::vector<EdgeData>                                          &edges,
    const std::unordered_map<CellId, std::vector<size_t>, CellIdHash>    &cell_edges,
    const std::unordered_map<size_t, std::vector<int64_t>>              &frozen_boundaries,
    const MagmaLattice                                                   &lattice,
    const MicronTables                                                   &um,
    int64_t horizon_um, int64_t min_h_um,
    int &out_windows, int &out_terms, int64_t &out_penalty_ub, int64_t &out_pair_cost)
{
    ortsat::LinearExpr penalty;
    out_windows = out_terms = 0;
    out_penalty_ub = out_pair_cost = 0;

    // ---- the layer grid for this block ----
    std::vector<int64_t> grid;
    grid.push_back(um.bottom_um[block.z_start_layer]);
    for (int L = block.z_start_layer; L <= block.z_end_layer; ++L)
        grid.push_back(um.top_um[L]);
    const int G = int(grid.size());
    if (G < 2)
        return penalty;

    std::unordered_map<int64_t, int> grid_index;
    for (int i = 0; i < G; ++i)
        grid_index[grid[i]] = i;

    // Window [i, win_end[i]) spans as close to horizon_um as the real layer heights allow, so
    // variable layer height is handled exactly rather than assumed away.
    std::vector<int> win_end(G);
    for (int i = 0, j = 0; i < G; ++i) {
        if (j < i) j = i;
        while (j < G && grid[j] - grid[i] < horizon_um) ++j;
        win_end[i] = j;
    }

    // Cost of one fully-coincident pair: it is charged once per window spanning its position,
    // so the widest window in grid steps. Measured in steps, not microns, so it stays correct
    // under variable layer height. Used to size the objective's gap allowance.
    for (int i = 0; i < G; ++i)
        out_pair_cost = std::max<int64_t>(out_pair_cost, win_end[i] - i);

    // ---- one-hot channelling: which grid position each boundary occupies ----
    // Sum(pos) == active, so an inactive tube contributes nothing to any count and no
    // enforcement literal is needed anywhere in the penalty path below.
    struct Flags { int lo, hi; std::vector<ortsat::BoolVar> pos; };   // grid range [lo, hi)
    std::vector<std::array<Flags, 2>> flags(segs.size());             // [slot][0=start,1=end]

    for (size_t si = 0; si < segs.size(); ++si) {
        for (int w = 0; w < 2; ++w) {
            const SlotRange r = boundary_range(segs[si], w == 1, min_h_um);
            int lo = G, hi = 0;
            for (int i = 0; i < G; ++i) {
                if (grid[i] < r.lo || grid[i] > r.hi) continue;
                lo = std::min(lo, i);
                hi = std::max(hi, i + 1);
            }
            if (lo >= hi) { flags[si][w].lo = flags[si][w].hi = 0; continue; }

            Flags &f = flags[si][w];
            f.lo = lo; f.hi = hi;
            f.pos.reserve(hi - lo);
            ortsat::LinearExpr sum, weighted;
            for (int i = lo; i < hi; ++i) {
                ortsat::BoolVar b = model.NewBoolVar();
                f.pos.push_back(b);
                sum      += b;
                weighted += ortsat::LinearExpr::Term(b, grid[i]);
            }
            model.AddEquality(sum, segs[si].active);
            model.AddEquality(w == 1 ? segs[si].end : segs[si].start, weighted)
                 .OnlyEnforceIf(segs[si].active);
        }
    }

    // ---- per cell: who feeds this cell's counter ----
    // Its own tubes and its neighbours'. A cell's own abutment seams are counted too; the
    // tube-count term already prices them.
    for (const CellId &cell : block.cells) {
        std::unordered_set<size_t> nbhd_edges;
        auto add_edges = [&](const CellId &c) {
            auto it = cell_edges.find(c);
            if (it == cell_edges.end()) return;
            for (size_t ei : it->second) nbhd_edges.insert(ei);
        };
        add_edges(cell);
        for (const CellId &nbr : lattice.neighbors(cell))
            add_edges(nbr);

        // Arrivals per grid position: decision boundaries as literals, frozen ones as counts.
        std::vector<ortsat::LinearExpr> arrivals(G);
        std::vector<int>                arrival_ub(G, 0);
        bool any = false;

        for (size_t si = 0; si < segs.size(); ++si) {
            if (!nbhd_edges.count(segs[si].edge_idx)) continue;
            for (int w = 0; w < 2; ++w) {
                const Flags &f = flags[si][w];
                for (int i = f.lo; i < f.hi; ++i) {
                    arrivals[i] += f.pos[i - f.lo];
                    ++arrival_ub[i];
                    any = true;
                }
            }
        }
        // Frozen boundaries over the same edge set, so each counts once, whichever cell holds it.
        for (size_t ei : nbhd_edges) {
            auto it = frozen_boundaries.find(ei);
            if (it == frozen_boundaries.end()) continue;
            for (int64_t z : it->second) {
                auto g = grid_index.find(z);
                if (g == grid_index.end()) continue;
                arrivals[g->second] += 1;
                ++arrival_ub[g->second];
                any = true;
            }
        }

        if (!any) continue;

        // Prefix sums, so a window count is one subtraction rather than a sum over its width.
        int total_ub = 0;
        for (int i = 0; i < G; ++i) total_ub += arrival_ub[i];

        std::vector<ortsat::IntVar> prefix(G);
        ortsat::LinearExpr running;
        for (int i = 0; i < G; ++i) {
            prefix[i] = model.NewIntVar(operations_research::Domain(0, total_ub));
            running += arrivals[i];
            model.AddEquality(prefix[i], running);
        }

        // One window per grid position, hinge-stacked.
        for (int i = 0; i < G; ++i) {
            const int j = win_end[i];
            if (j <= i + 1) continue;                 // a window holding one position cannot crowd

            int win_ub = 0;
            for (int k = i; k < j; ++k) win_ub += arrival_ub[k];
            if (win_ub < 2) continue;                 // nothing here could ever collide

            // count = prefix[j-1] - prefix[i-1]
            ortsat::LinearExpr count = prefix[j - 1];
            if (i > 0) count -= prefix[i - 1];

            ++out_windows;
            for (int lvl = 1; lvl <= STAGGER_HINGE_LEVELS && lvl < win_ub; ++lvl) {
                const int64_t ub = win_ub - lvl;
                ortsat::IntVar pen = model.NewIntVar(operations_research::Domain(0, ub));
                model.AddGreaterOrEqual(pen, count - lvl);
                penalty += pen;
                out_penalty_ub += ub;
                ++out_terms;
            }
        }
    }

    return penalty;
}

} // namespace

// ============================================================================
// solve_block — build and solve CP-SAT model for one block
// ============================================================================

BlockResult MagmaTubeSolver::solve_block(const Block &block) const
{
    using namespace operations_research::sat;

    const int64_t min_h_um = llround(m_min_h_mm * 1000.0);
    const int64_t max_h_um = llround(m_max_h_mm * 1000.0);
    const int64_t horizon_um = stagger_horizon_um(min_h_um);

    CpModelBuilder model;

    // Track intervals per cell for NoOverlap
    std::unordered_map<CellId, std::vector<IntervalVar>, CellIdHash> cell_intervals;
    // Per-cell committed coverage, accumulated only over runs this model actually represents.
    std::unordered_map<CellId, int64_t, CellIdHash> cell_floor;

    std::vector<SegVars> all_segments;
    // Summed per run below; becomes a hard cap on coverage_expr once the model is built.
    int64_t coverage_ceiling_um = 0;

    // ------------------------------------------------------------------
    // 1. Create decision variables for each edge's runs
    //
    // All micron-space variables use discrete domains derived from a
    // unified layer boundary list. Since bottom_um[L+1] == top_um[L],
    // each boundary is both the end of one layer and the start of the
    // next. One sorted list drives all domains:
    //   position (start/end) ∈ boundaries
    //   size    ∈ {b[j] - b[i] | j > i} ∩ [min_h, max_h]
    //   contrib ∈ {0} ∪ sizes
    // ------------------------------------------------------------------
    for (size_t ei : block.edge_indices) {
        const EdgeData &ed = m_edges[ei];

        for (const Run &run : ed.runs) {
            // Skip runs that don't overlap block Z range
            if (run.end_um <= block.z_start_um || run.start_um >= block.z_end_um)
                continue;

            // Restrict to the intersection of run and block Z range.
            int eff_start = std::max(run.start_layer, block.z_start_layer);
            int eff_end   = std::min(run.end_layer, block.z_end_layer);
            if (eff_start > eff_end)
                continue;

            int64_t eff_h = m_um.top_um[eff_end] - m_um.bottom_um[eff_start];
            if (eff_h < min_h_um)
                continue; // Clipped range too short for a tube

            const int64_t run_start_um = m_um.bottom_um[eff_start];
            const int64_t run_end_um   = m_um.top_um[eff_end];

            // Slot count: how many tubes fit in the run, independent of what is committed.
            const int K = std::max<int>(2, int((run_end_um - run_start_um) / min_h_um));

            // Unified layer boundary list for this run
            std::vector<int64_t> boundaries;
            boundaries.push_back(m_um.bottom_um[eff_start]);
            for (int L = eff_start; L <= eff_end; ++L)
                boundaries.push_back(m_um.top_um[L]);
            // boundaries is sorted: each value is the end of one layer
            // and the start of the next

            // Feasible sizes: all boundary-pair differences in height range
            std::set<int64_t> size_set;
            for (size_t i = 0; i < boundaries.size(); ++i)
                for (size_t j = i + 1; j < boundaries.size(); ++j) {
                    int64_t s = boundaries[j] - boundaries[i];
                    if (s > max_h_um) break; // sorted, no point continuing
                    if (s >= min_h_um)
                        size_set.insert(s);
                }

            if (size_set.empty())
                continue; // No valid tube height for this run

            // The per-cell coverage floor counts only committed segments inside a modelled
            // run: anything else has no variable to satisfy it and would make the block
            // infeasible.
            for (const CommittedSegment &cs : m_committed[ei])
                if (cs.start_um >= run_start_um && cs.end_um <= run_end_um) {
                    const int64_t h = cs.end_um - cs.start_um;
                    cell_floor[ed.edge.a] += h;
                    cell_floor[ed.edge.b] += h;
                }

            // Upper bound on this run's coverage: at most L/min_size disjoint tubes of at most
            // max_size each, capped by L. Tighter than L when no tube count tiles the run
            // (e.g. min_h 2.0, max_h 3.4, L 3.7). Never excludes a legal solution.
            {
                const int64_t L        = run_end_um - run_start_um;
                const int64_t min_size = *size_set.begin();
                const int64_t max_size = *size_set.rbegin();
                coverage_ceiling_um += std::min(L, (L / min_size) * max_size);
            }

            std::vector<int64_t> size_vals(size_set.begin(), size_set.end());
            std::vector<int64_t> contrib_vals = {0};
            contrib_vals.insert(contrib_vals.end(), size_vals.begin(), size_vals.end());

            auto boundary_dom = operations_research::Domain::FromValues(boundaries);
            auto size_dom     = operations_research::Domain::FromValues(size_vals);
            auto contrib_dom  = operations_research::Domain::FromValues(contrib_vals);

            BoolVar prev_active;
            bool has_prev = false;

            for (int k = 0; k < K; ++k) {
                BoolVar active = model.NewBoolVar();
                IntVar  start  = model.NewIntVar(boundary_dom);
                IntVar  end    = model.NewIntVar(boundary_dom);
                IntVar  size   = model.NewIntVar(size_dom);

                // size = end - start
                model.AddEquality(size, end - start);

                // Optional interval (half-open: [start, start+size) = [start, end))
                IntervalVar interval = model.NewOptionalIntervalVar(
                    start, size, end, active);

                // Add to BOTH cells' NoOverlap
                cell_intervals[ed.edge.a].push_back(interval);
                cell_intervals[ed.edge.b].push_back(interval);

                // Symmetry breaking within run
                if (has_prev) {
                    // Segment k can only be active if k-1 is active
                    model.AddImplication(active, prev_active);
                    // When both active, k starts at or after k-1's end. Must stay
                    // non-strict: positions are layer boundaries, so strict
                    // less-than would force a layer-high gap between abutting
                    // tubes and reject greedy's hint, which abuts them.
                    model.AddLessOrEqual(all_segments.back().end, start)
                        .OnlyEnforceIf({prev_active, active});
                }

                int64_t rs = boundaries.front();
                int64_t re = boundaries.back();
                all_segments.push_back({active, start, end, size,
                                        ei, rs, re, k, contrib_dom});
                prev_active = active;
                has_prev = true;
            }
        }
    }

    // ------------------------------------------------------------------
    // 2. Add frozen intervals from committed tubes outside the block
    // ------------------------------------------------------------------

    // Frozen boundaries for stagger penalty
    // Keyed by edge, so the stagger penalty can count each boundary once per neighbourhood.
    std::unordered_map<size_t, std::vector<int64_t>> frozen_boundaries;

    // Scan boundary edges via cell reverse lookup (avoids scanning all edges)
    std::unordered_set<size_t> seen_boundary_edges;
    for (const CellId &cell : block.cells) {
        auto it = m_cell_edges.find(cell);
        if (it == m_cell_edges.end()) continue;
        for (size_t ei : it->second) {
            if (!seen_boundary_edges.insert(ei).second) continue;
            const EdgeData &ed = m_edges[ei];
            bool a_in = block.cells.count(ed.edge.a) > 0;
            bool b_in = block.cells.count(ed.edge.b) > 0;

            // Skip edges fully inside (decision edges) or fully outside (irrelevant)
            if (a_in == b_in)
                continue;

            const CellId &inside_cell = a_in ? ed.edge.a : ed.edge.b;

            for (const CommittedSegment &seg : m_committed[ei]) {
                // Clip to block Z range
                int64_t clip_start = std::max(seg.start_um, block.z_start_um);
                int64_t clip_end   = std::min(seg.end_um, block.z_end_um);
                if (clip_start >= clip_end)
                    continue;

                IntervalVar frozen = model.NewFixedSizeIntervalVar(
                    clip_start, clip_end - clip_start);
                cell_intervals[inside_cell].push_back(frozen);

                // Record boundaries for stagger (use actual, not clipped)
                frozen_boundaries[ei].push_back(seg.start_um);
                frozen_boundaries[ei].push_back(seg.end_um);
            }
        }
    }

    // Decision-edge segments that extend outside the block's Z range are frozen; those fully
    // inside become warm-start hints and are re-optimized.
    for (size_t ei : block.edge_indices) {
        for (const CommittedSegment &seg : m_committed[ei]) {
            if (seg.start_um >= block.z_start_um && seg.end_um <= block.z_end_um)
                continue; // Fully inside — warm start, not frozen

            int64_t clip_start = std::max(seg.start_um, block.z_start_um);
            int64_t clip_end   = std::min(seg.end_um, block.z_end_um);
            if (clip_start >= clip_end)
                continue;

            IntervalVar frozen = model.NewFixedSizeIntervalVar(
                clip_start, clip_end - clip_start);
            const EdgeData &ed = m_edges[ei];
            cell_intervals[ed.edge.a].push_back(frozen);
            cell_intervals[ed.edge.b].push_back(frozen);

            frozen_boundaries[ei].push_back(seg.start_um);
            frozen_boundaries[ei].push_back(seg.end_um);
        }
    }

    // No decision variables (e.g. every run clipped too short by the Z range): nothing to solve.
    if (all_segments.empty())
        return BlockResult{BlockOutcome::NoWork, {}};

    // ------------------------------------------------------------------
    // 3. NoOverlap per cell
    // ------------------------------------------------------------------
    for (auto &[cell, intervals] : cell_intervals) {
        if (intervals.size() > 1)
            model.AddNoOverlap(intervals);
    }

    // ------------------------------------------------------------------
    // 4. Contribution variables
    //
    // A slot's contribution is its size when active and 0 when not. The objective maximises
    // their sum and pays for active_expr; the per-cell floors below constrain their sums per
    // cell.
    // ------------------------------------------------------------------

    LinearExpr coverage_expr;  // Σ contrib, in µm
    LinearExpr active_expr;    // Σ active
    std::vector<IntVar> seg_contrib;
    seg_contrib.reserve(all_segments.size());

    for (const auto &seg : all_segments) {
        IntVar contrib = model.NewIntVar(seg.contrib_dom);
        model.AddEquality(contrib, seg.size).OnlyEnforceIf(seg.active);
        model.AddEquality(contrib, 0).OnlyEnforceIf(seg.active.Not());
        coverage_expr += contrib;
        active_expr   += seg.active;
        seg_contrib.push_back(contrib);
    }

    // ------------------------------------------------------------------
    // 4b. Per-cell coverage floors
    // ------------------------------------------------------------------
    //
    // A cell's committed coverage is a floor the solve cannot go below. The objective weight
    // only guarantees coverage for a solve that reaches optimality; the floor also holds when
    // a block hits its time limit. Per cell, because an aggregate floor would allow emptying
    // one cell to over-fill its neighbour.
    {
        std::unordered_map<CellId, ortsat::LinearExpr, CellIdHash> cell_cov;
        for (size_t si = 0; si < all_segments.size(); ++si) {
            const CellEdge &e = m_edges[all_segments[si].edge_idx].edge;
            cell_cov[e.a] += seg_contrib[si];
            cell_cov[e.b] += seg_contrib[si];
        }
        for (const auto &cf : cell_floor) {
            auto it = cell_cov.find(cf.first);
            if (it != cell_cov.end() && cf.second > 0)
                model.AddGreaterOrEqual(it->second, cf.second);
        }
    }

    // ------------------------------------------------------------------
    // 4c. Objective: coverage first, then tube count + stagger penalty
    // ------------------------------------------------------------------
    //
    // A single solve, so stagger chooses among the coverage-optimal solutions rather than
    // inheriting an arbitrary one. Coverage is in integer microns, and its weight exceeds the
    // largest possible sum of the secondary terms (slot count + stagger_ub), so one micron of
    // fill outscores any secondary trade. Tube count and stagger share the secondary tier at
    // equal unit weight.
    int     stagger_windows = 0, stagger_terms = 0;
    int64_t stagger_ub = 0, stagger_pair_cost = 0;
    LinearExpr stagger;

    if (m_mode == MagmaTubeSolverMode::CoverageStagger && horizon_um > 0)
        stagger = build_stagger_penalty(
            model, block, all_segments, m_edges, m_cell_edges, frozen_boundaries,
            m_lattice, m_um, horizon_um, min_h_um,
            stagger_windows, stagger_terms, stagger_ub, stagger_pair_cost);

    // Redundant cuts linking coverage to tube count, which the LP relaxation cannot derive
    // cheaply from the per-segment relations; they shorten the optimality proof.
    model.AddGreaterOrEqual(coverage_expr, min_h_um * active_expr);
    model.AddLessOrEqual(coverage_expr, max_h_um * active_expr);

    // Constant ceiling on coverage. The relaxation escapes the cuts above by raising
    // fractional activity (NoOverlap gives no linear cut), so without this the dual bound
    // overstates coverage and the gap never closes to stagger size.
    if (coverage_ceiling_um > 0)
        model.AddLessOrEqual(coverage_expr, coverage_ceiling_um);

    // W is also inflated by the absolute gap limit set below, so a micron of coverage is worth
    // more than the gap and stopping early can never cost coverage; the secondary terms land
    // within STAGGER_GAP_PAIRS coincident pairs of optimal. stagger_ub is the exact sum of the
    // hinge variables' upper bounds.
    const int64_t secondary_max = int64_t(all_segments.size()) + stagger_ub;
    const int64_t stagger_gap   = stagger_terms > 0 ? STAGGER_GAP_PAIRS * stagger_pair_cost : 0;
    const int64_t W_COVERAGE    = secondary_max + stagger_gap + 1;

    LinearExpr objective = W_COVERAGE * coverage_expr - active_expr;
    if (stagger_terms > 0)
        objective -= stagger;

    model.Maximize(objective);

    // ------------------------------------------------------------------
    // 5. Complete solution hint (warm start)
    // ------------------------------------------------------------------
    {
        std::vector<bool> hinted(all_segments.size(), false);

        // Collect and sort committed segments per edge (within Z range)
        std::unordered_map<size_t, std::vector<const CommittedSegment *>> edge_hints;
        for (size_t ei : block.edge_indices) {
            for (const CommittedSegment &seg : m_committed[ei]) {
                if (seg.start_um >= block.z_start_um && seg.end_um <= block.z_end_um)
                    edge_hints[ei].push_back(&seg);
            }
        }
        for (auto &[ei, hints] : edge_hints)
            std::sort(hints.begin(), hints.end(),
                      [](const CommittedSegment *a, const CommittedSegment *b) {
                          return a->start_um < b->start_um;
                      });

        // Match committed segments, in order, to the leading slots of the run containing
        // them, which keeps hinted activation a prefix within each run.
        for (auto &[ei, hints] : edge_hints) {
            if (hints.empty()) continue;
            size_t ci = 0;
            for (size_t si = 0; si < all_segments.size(); ++si) {
                if (ci >= hints.size()) break;
                const auto &sv = all_segments[si];
                if (sv.edge_idx != ei) continue;
                const CommittedSegment *cs = hints[ci];
                if (cs->start_um >= sv.run_start_um && cs->end_um <= sv.run_end_um) {
                    model.AddHint(sv.active, true);
                    model.AddHint(sv.start, cs->start_um);
                    model.AddHint(sv.end, cs->end_um);
                    hinted[si] = true;
                    ++ci;
                }
            }
        }

        // Hint the remaining slots inactive. CP-SAT silently discards the whole hint if any
        // part of it violates a hard constraint, so hints must stay feasible.
        for (size_t si = 0; si < all_segments.size(); ++si) {
            if (!hinted[si])
                model.AddHint(all_segments[si].active, false);
        }
    }

    // ------------------------------------------------------------------
    // 6. Solve
    // ------------------------------------------------------------------
    SatParameters params;
    params.set_max_time_in_seconds(m_per_block_timeout);
    params.set_num_workers(m_cpsat_workers);
    // Stop once no solution can be more than stagger_gap better; W_COVERAGE covers that
    // allowance, so this never costs coverage. Must be absolute: a relative gap on a
    // coverage-dominated objective would admit losing whole tubes.
    if (stagger_gap > 0)
        params.set_absolute_gap_limit(double(stagger_gap));
    // Default linearization_level performs better with these discrete domains.
    // Do not set repair_hint or hint_conflict_limit: they trigger a CP-SAT fixed_search
    // crash (integer_search.cc:1217).

    auto proto = model.Build();
    std::string validation_error = ValidateCpModel(proto);
    if (!validation_error.empty()) {
        BOOST_LOG_TRIVIAL(error) << "MagmaTubeSolver: MODEL INVALID: "
            << validation_error.substr(0, 200);
    }

    operations_research::sat::Model sat_model;
    sat_model.Add(NewSatParameters(params));

    // Cancellation is checked between blocks only. A mid-solve abort via
    // TimeLimit::RegisterExternalBooleanAsLimit conflicts with NewSatParameters.

    CpSolverResponse response = SolveCpModel(proto, &sat_model);

    // ------------------------------------------------------------------
    // 7. Extract results
    // ------------------------------------------------------------------
    BlockResult result;
    if (response.status() == CpSolverStatus::OPTIMAL ||
        response.status() == CpSolverStatus::FEASIBLE) {
        result.outcome = BlockOutcome::Solved;
        for (const auto &seg : all_segments) {
            if (SolutionBooleanValue(response, seg.active)) {
                int64_t s = SolutionIntegerValue(response, seg.start);
                int64_t e = SolutionIntegerValue(response, seg.end);
                result.segments.push_back({seg.edge_idx, {s, e}});
            }
        }
        // Per-block coverage comparison
        int64_t solver_cov = 0;
        for (const auto &[eidx, seg] : result.segments)
            solver_cov += seg.end_um - seg.start_um;
        int64_t committed_cov = 0;
        for (size_t ei : block.edge_indices)
            for (const CommittedSegment &cs : m_committed[ei])
                if (cs.start_um >= block.z_start_um && cs.end_um <= block.z_end_um)
                    committed_cov += cs.end_um - cs.start_um;

        // Diagnostic only: the objective is coverage-dominated, so this says little about stagger.
        const double gap_pct =
            (response.best_objective_bound() > 0)
                ? (1.0 - response.objective_value() / response.best_objective_bound()) * 100.0
                : 0.0;

        BOOST_LOG_TRIVIAL(debug) << "MagmaTubeSolver: block solved "
            << (response.status() == CpSolverStatus::OPTIMAL ? "OPTIMAL" : "FEASIBLE")
            << ", obj=" << response.objective_value()
            << ", gap=" << gap_pct << "%"
            << ", segs=" << result.segments.size()
            << "/" << all_segments.size()
            << ", wall=" << std::fixed << std::setprecision(2)
            << response.wall_time() << "s"
            << " | cov=" << solver_cov/1000.0
            << "mm committed=" << committed_cov/1000.0 << "mm"
            << " ceil=" << coverage_ceiling_um/1000.0 << "mm"
            << " | windows=" << stagger_windows << " terms=" << stagger_terms
            << " " << response.solution_info();

        // The per-cell floors should make this unreachable; checked rather than trusted.
        if (solver_cov < committed_cov)
            BOOST_LOG_TRIVIAL(error) << "  COVERAGE DROP: lost "
                << (committed_cov - solver_cov)/1000.0 << "mm against what was committed";
    } else {
        // Interval counts per cell, for diagnosing an unsolvable block.
        int total_intervals = 0;
        int max_intervals_per_cell = 0;
        for (const auto &[cell, intervals] : cell_intervals) {
            int n = static_cast<int>(intervals.size());
            if (n > max_intervals_per_cell) max_intervals_per_cell = n;
            total_intervals += n;
        }
        BOOST_LOG_TRIVIAL(warning) << "MagmaTubeSolver: block solve status="
            << static_cast<int>(response.status())
            << " (" << (response.status() == CpSolverStatus::INFEASIBLE ? "INFEASIBLE" :
                        response.status() == CpSolverStatus::MODEL_INVALID ? "MODEL_INVALID" :
                        response.status() == CpSolverStatus::UNKNOWN ? "UNKNOWN" : "OTHER")
            << "), edges=" << block.edge_indices.size()
            << ", segments=" << all_segments.size()
            << ", cells=" << block.cells.size()
            << ", cell_interval_lists=" << cell_intervals.size()
            << ", total_intervals=" << total_intervals
            << ", max_per_cell=" << max_intervals_per_cell;
    }

    return result;
}

// ============================================================================
// solve_pass
// ============================================================================

void MagmaTubeSolver::solve_pass(const std::vector<Block> &blocks,
                                 std::function<void(int)> block_done_fn)
{
    auto t_start = std::chrono::high_resolution_clock::now();

    // Sequential: each block commits before the next, so overlapping blocks see the latest
    // state. CP-SAT parallelises within a block.
    for (size_t i = 0; i < blocks.size(); ++i) {
        // Re-divide the remaining budget so time unused by easy blocks passes to hard ones.
        // The floor can push the total past the budget.
        m_per_block_timeout = std::max(0.5, m_budget_left_sec / std::max(1, m_blocks_left));

        const auto t_block = std::chrono::high_resolution_clock::now();
        commit_result(blocks[i], solve_block(blocks[i]));
        const double spent = std::chrono::duration<double>(
            std::chrono::high_resolution_clock::now() - t_block).count();

        m_budget_left_sec = std::max(0.0, m_budget_left_sec - spent);
        if (m_blocks_left > 0) --m_blocks_left;

        if (block_done_fn) block_done_fn(static_cast<int>(i));
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

    BOOST_LOG_TRIVIAL(debug) << "MagmaTubeSolver: pass: "
        << blocks.size() << " blocks, " << ms << "ms";
}

// ============================================================================
// commit_result
// ============================================================================

void MagmaTubeSolver::commit_result(const Block &block, const BlockResult &result)
{
    if (result.outcome == BlockOutcome::NoWork)
        return;
    if (result.outcome == BlockOutcome::Failed) {
        ++m_failed_blocks;
        return;
    }

    // Replace only segments fully inside this block's Z range: a block modifies only what it
    // can fully see. Segments extending outside were frozen in its model.
    for (size_t ei : block.edge_indices) {
        auto &segs = m_committed[ei];
        segs.erase(std::remove_if(segs.begin(), segs.end(),
            [&](const CommittedSegment &s) {
                return s.start_um >= block.z_start_um &&
                       s.end_um   <= block.z_end_um;
            }), segs.end());
    }
    for (const auto &[edge_idx, seg] : result.segments)
        m_committed[edge_idx].push_back(seg);
}

// ============================================================================
// extract_results — convert micron-space assignments to UTubePair format
// ============================================================================

void MagmaTubeSolver::extract_results(
    std::vector<UTubePair> &out_pairs,
    std::unordered_map<CellId, std::vector<int>, CellIdHash> &out_index) const
{
    out_pairs.clear();
    out_index.clear();

    for (size_t ei = 0; ei < m_edges.size(); ++ei) {
        const EdgeData &ed = m_edges[ei];
        for (const CommittedSegment &seg : m_committed[ei]) {
            auto start_it = m_um.bottom_to_layer.find(seg.start_um);
            auto end_it   = m_um.top_to_layer.find(seg.end_um);
            if (start_it == m_um.bottom_to_layer.end() ||
                end_it   == m_um.top_to_layer.end()) {
                BOOST_LOG_TRIVIAL(warning) << "MagmaTubeSolver: invalid micron boundary "
                    << seg.start_um << "-" << seg.end_um;
                continue;
            }

            UTubePair pair;
            pair.cell_a           = ed.edge.a;
            pair.cell_b           = ed.edge.b;
            pair.pair_start_layer = start_it->second;
            pair.pair_end_layer   = end_it->second;
            pair.volume_mm3       = 0.0; // computed later by MagmaTubeMap
            pair.window_end_z     = 0.0; // computed later by MagmaTubeMap

            int pair_idx = static_cast<int>(out_pairs.size());
            out_pairs.push_back(pair);
            out_index[ed.edge.a].push_back(pair_idx);
            out_index[ed.edge.b].push_back(pair_idx);
        }
    }

    // Ensure all cells have an entry (empty = solid fill)
    for (const auto &[cell, _] : m_cells) {
        if (out_index.find(cell) == out_index.end())
            out_index[cell]; // insert empty vector
    }

    BOOST_LOG_TRIVIAL(info) << "MagmaTubeSolver: extracted " << out_pairs.size() << " pairs";
}

// ============================================================================
// validate_and_report — per-label validation + diagnostics
// ============================================================================

void MagmaTubeSolver::validate_and_report(const char *label) const
{
    ValidationResult result;

    // Build per-cell segment list for overlap checking
    std::unordered_map<CellId, std::vector<std::pair<int, int>>,
                       CellIdHash> cell_segments; // (start_layer, end_layer)

    int total_segs = 0;
    for (size_t ei = 0; ei < m_edges.size(); ++ei) {
        const EdgeData &ed = m_edges[ei];
        for (const CommittedSegment &seg : m_committed[ei]) {
            ++total_segs;
            auto start_it = m_um.bottom_to_layer.find(seg.start_um);
            auto end_it   = m_um.top_to_layer.find(seg.end_um);
            if (start_it == m_um.bottom_to_layer.end() ||
                end_it   == m_um.top_to_layer.end()) {
                if (++result.bad_range <= 3)
                    BOOST_LOG_TRIVIAL(warning) << label << ": invalid micron boundary "
                        << seg.start_um << "-" << seg.end_um;
                continue;
            }

            int sL = start_it->second;
            int eL = end_it->second;

            // 1. Valid layer range
            if (sL < 0 || eL >= m_num_layers || sL > eL) {
                if (++result.bad_range <= 3)
                    BOOST_LOG_TRIVIAL(warning) << label << ": tube invalid range L"
                        << sL << "-" << eL;
                continue;
            }

            // 2. Height bounds
            double h_mm = m_layer_data[eL].print_z
                        - (m_layer_data[sL].print_z - m_layer_data[sL].height);
            if (h_mm < m_min_h_mm - 0.01) {
                if (++result.bad_short <= 3)
                    BOOST_LOG_TRIVIAL(warning) << label << ": tube too short: "
                        << h_mm << "mm (min=" << m_min_h_mm << ") L" << sL << "-" << eL;
            }
            if (h_mm > m_max_h_mm + 0.01) {
                if (++result.bad_long <= 3)
                    BOOST_LOG_TRIVIAL(warning) << label << ": tube too long: "
                        << h_mm << "mm (max=" << m_max_h_mm << ") L" << sL << "-" << eL;
            }

            // 3. Both cells present at every layer
            auto it_a = m_cells.find(ed.edge.a);
            auto it_b = m_cells.find(ed.edge.b);
            if (it_a != m_cells.end() && it_b != m_cells.end()) {
                for (int L = sL; L <= eL; ++L) {
                    if (!it_a->second.present(L) || !it_b->second.present(L)) {
                        if (++result.bad_presence <= 3)
                            BOOST_LOG_TRIVIAL(warning) << label << ": cell not present at L"
                                << L << " (range L" << sL << "-" << eL << ")";
                        break;
                    }
                }
            }

            cell_segments[ed.edge.a].push_back({sL, eL});
            cell_segments[ed.edge.b].push_back({sL, eL});
        }
    }

    // 4. Per-cell overlap
    for (const auto &[cell, segs] : cell_segments) {
        for (size_t i = 0; i < segs.size(); ++i) {
            for (size_t j = i + 1; j < segs.size(); ++j) {
                if (segs[i].first <= segs[j].second &&
                    segs[j].first <= segs[i].second) {
                    ++result.overlap_cell;
                    if (result.overlap_cell <= 3)
                        BOOST_LOG_TRIVIAL(warning) << label << ": cell ("
                            << cell.a << "," << cell.b << "," << cell.c
                            << ") overlap [L" << segs[i].first << "-" << segs[i].second
                            << "] vs [L" << segs[j].first << "-" << segs[j].second << "]";
                }
            }
        }
    }

    // 5. Per-edge overlap
    for (size_t ei = 0; ei < m_edges.size(); ++ei) {
        const auto &segs = m_committed[ei];
        for (size_t i = 0; i < segs.size(); ++i) {
            for (size_t j = i + 1; j < segs.size(); ++j) {
                if (segs[i].start_um < segs[j].end_um &&
                    segs[j].start_um < segs[i].end_um) {
                    ++result.overlap_edge;
                    if (result.overlap_edge <= 3)
                        BOOST_LOG_TRIVIAL(warning) << label << ": edge " << ei
                            << " overlap [" << segs[i].start_um << "-" << segs[i].end_um
                            << "] vs [" << segs[j].start_um << "-" << segs[j].end_um << "]";
                }
            }
        }
    }

    // Summary
    if (result.total_issues())
        BOOST_LOG_TRIVIAL(warning) << label << ": " << result.total_issues() << " issues ("
            << result.bad_short << " short, " << result.bad_long << " long, "
            << result.bad_range << " bad-range, "
            << result.bad_presence << " not-present, "
            << result.overlap_cell << " cell-overlap, "
            << result.overlap_edge << " edge-overlap)";
    else
        BOOST_LOG_TRIVIAL(info) << label << ": all " << total_segs << " segments OK";

    report_coverage(cell_segments, label);
    report_lengths(label);
    report_stagger(label);
}

// ============================================================================
// report_coverage — how much of each cell's presence range ended up in a tube
// ============================================================================

void MagmaTubeSolver::report_coverage(
    const std::unordered_map<CellId, std::vector<std::pair<int, int>>,
                             CellIdHash> &cell_segments,
    const char *label) const
{
    double total_presence_um = 0, total_covered_um = 0;
    int cells_below_50 = 0, cells_below_25 = 0, cells_zero = 0;
    for (const auto &[cell, presence] : m_cells) {
        // Only layers the cell is present on: presence can have gaps (see CellPresence),
        // which no tube may cross.
        int64_t presence_um = 0;
        for (int L = presence.first_layer; L <= presence.last_layer; ++L)
            if (presence.present(L))
                presence_um += m_um.top_um[L] - m_um.bottom_um[L];
        if (presence_um <= 0) continue;

        int64_t covered_um = 0;
        auto it = cell_segments.find(cell);
        if (it != cell_segments.end()) {
            // A tube counts for both of its cells, since it reinforces both.
            for (const auto &[sL, eL] : it->second)
                covered_um += m_um.top_um[eL] - m_um.bottom_um[sL];
        }
        double pct = double(covered_um) / double(presence_um);
        total_presence_um += presence_um;
        total_covered_um  += covered_um;
        // Cumulative bands: a cell under 25% is also counted under 50%.
        if (covered_um == 0) ++cells_zero;
        if (pct < 0.25)      ++cells_below_25;
        if (pct < 0.50)      ++cells_below_50;
    }
    double overall_pct = total_presence_um > 0
        ? total_covered_um / total_presence_um * 100.0 : 0.0;
    BOOST_LOG_TRIVIAL(info) << label << " COVERAGE: " << std::fixed << std::setprecision(1)
        << overall_pct << "% overall, "
        << cells_zero << " cells unfilled, "
        << cells_below_25 << " cells <25%, "
        << cells_below_50 << " cells <50%";
}

// ============================================================================
// report_lengths — tube-height distribution
// ============================================================================
// Uniform tube heights put every boundary on one regular grid, so the height spread
// predicts the stagger numbers.

void MagmaTubeSolver::report_lengths(const char *label) const
{
    std::vector<int64_t> heights;
    for (const auto &segs : m_committed)
        for (const CommittedSegment &seg : segs)
            heights.push_back(seg.end_um - seg.start_um);

    if (heights.empty()) {
        BOOST_LOG_TRIVIAL(info) << label << " LENGTHS: no tubes";
        return;
    }

    std::sort(heights.begin(), heights.end());
    const size_t  n      = heights.size();
    const int64_t min_um = heights.front();
    const int64_t max_um = heights.back();
    const int64_t med_um = heights[n / 2];
    int64_t sum_um = 0;
    for (int64_t h : heights) sum_um += h;

    // Distinct heights, most common first. Tubes snap to layer boundaries, so the list is short.
    std::vector<std::pair<int64_t, size_t>> hist; // (height_um, count)
    for (int64_t h : heights) {
        if (!hist.empty() && hist.back().first == h) ++hist.back().second;
        else hist.push_back({h, 1});
    }
    const size_t distinct = hist.size();
    std::sort(hist.begin(), hist.end(),
              [](const auto &x, const auto &y) { return x.second > y.second; });

    std::ostringstream top;
    for (size_t i = 0; i < hist.size() && i < 5; ++i)
        top << (i ? ", " : "") << std::fixed << std::setprecision(2)
            << hist[i].first / 1000.0 << "mm x" << hist[i].second;

    // Same-edge abutments: consecutive tubes on one edge that meet exactly. Non-zero on the
    // GREEDY line but zero on CPSAT suggests the warm-start hint is being rejected (see the
    // slot ordering constraint in solve_block).
    size_t abutments = 0;
    std::vector<CommittedSegment> ordered;
    for (const auto &segs : m_committed) {
        if (segs.size() < 2) continue;
        ordered.assign(segs.begin(), segs.end());
        std::sort(ordered.begin(), ordered.end(),
                  [](const CommittedSegment &x, const CommittedSegment &y) {
                      return x.start_um < y.start_um;
                  });
        for (size_t i = 1; i < ordered.size(); ++i)
            if (ordered[i].start_um == ordered[i - 1].end_um) ++abutments;
    }

    BOOST_LOG_TRIVIAL(info) << label << " LENGTHS: " << n << " tubes, "
        << std::fixed << std::setprecision(2)
        << "min=" << min_um / 1000.0
        << " median=" << med_um / 1000.0
        << " mean=" << double(sum_um) / double(n) / 1000.0
        << " max=" << max_um / 1000.0 << "mm"
        << " (cap " << m_max_h_mm << "mm)"
        << " | " << distinct << " distinct: " << top.str()
        << (distinct > 5 ? ", ..." : "")
        << " | " << abutments << " same-edge abutments";
}

// ============================================================================
// report_stagger — how far apart neighbouring tube boundaries sit
// ============================================================================
// Boundaries at the same Z across touching cells form a plane the part can split along.
// For each cell, gathers the boundaries of its own tubes (Ring-0) and its neighbours'
// (Ring-1), the same neighbourhood the stagger objective uses, and scores pairs closer than
// the horizon, split into avoidable and forced.

void MagmaTubeSolver::report_stagger(const char *label) const
{
    const int64_t horizon_um = stagger_horizon_um(llround(m_min_h_mm * 1000.0));
    if (horizon_um <= 0) {
        BOOST_LOG_TRIVIAL(warning) << label
            << " STAGGER: no horizon -- min tube height is zero, which should be impossible";
        return;
    }
    const int64_t tight_um = std::max<int64_t>(1, horizon_um / 2);

    // Each segment contributes two boundaries to both of its cells. A boundary carries its
    // segment id so duplicates (one segment reached via several cells) are removed without
    // merging two different tubes that break at the same Z, and its edge so a pair can be
    // classified.
    struct Boundary {
        int64_t  z_um;
        uint32_t seg_uid;
        size_t   edge_idx;
    };
    std::unordered_map<CellId, std::vector<Boundary>, CellIdHash> cell_bounds;
    uint32_t next_uid = 0;
    for (size_t ei = 0; ei < m_edges.size(); ++ei) {
        for (const CommittedSegment &seg : m_committed[ei]) {
            const uint32_t uid = next_uid++;
            for (const CellId &c : {m_edges[ei].edge.a, m_edges[ei].edge.b}) {
                cell_bounds[c].push_back({seg.start_um, uid, ei});
                cell_bounds[c].push_back({seg.end_um,   uid, ei});
            }
        }
    }
    if (next_uid == 0) {
        BOOST_LOG_TRIVIAL(info) << label << " STAGGER: no tubes";
        return;
    }

    auto shares_cell = [&](size_t ei, size_t ej) {
        const CellEdge &A = m_edges[ei].edge, &B = m_edges[ej].edge;
        return A.a == B.a || A.a == B.b || A.b == B.a || A.b == B.b;
    };

    // Pairs are reported in two populations. Forced seams grow with coverage and no solver
    // choice removes them, so mixing them in would make the count track coverage, not stagger.
    //   AVOIDABLE -- the two tubes share no cell: independent columns breaking at the same Z.
    //   FORCED    -- the two tubes share a cell, whose column they tile, so they can only come
    //                close by abutting, which full coverage requires.
    // The stagger objective scores both populations.
    struct Bucket {
        size_t  pairs = 0, within_tight = 0, coincident = 0;
        int64_t shortfall_um = 0;
        int64_t min_gap_um = std::numeric_limits<int64_t>::max();
    };
    Bucket avoidable, forced;

    // A pair reachable from several cells' neighbourhoods is counted once.
    std::set<std::array<int64_t, 4>> seen_pairs;
    std::unordered_map<CellId, size_t, CellIdHash> cell_offences;

    std::vector<Boundary> nb;
    for (const auto &cell_entry : m_cells) {
        const CellId &cell = cell_entry.first;
        nb.clear();
        auto append = [&](const CellId &c) {
            auto it = cell_bounds.find(c);
            if (it != cell_bounds.end())
                nb.insert(nb.end(), it->second.begin(), it->second.end());
        };
        append(cell);
        for (const CellId &nbr : m_lattice.neighbors(cell))
            append(nbr);
        if (nb.size() < 2) continue;

        std::sort(nb.begin(), nb.end(), [](const Boundary &x, const Boundary &y) {
            return x.z_um != y.z_um ? x.z_um < y.z_um : x.seg_uid < y.seg_uid;
        });
        nb.erase(std::unique(nb.begin(), nb.end(),
                     [](const Boundary &x, const Boundary &y) {
                         return x.z_um == y.z_um && x.seg_uid == y.seg_uid;
                     }),
                 nb.end());

        // Sliding window: only boundaries within the horizon can score at all.
        for (size_t i = 0; i < nb.size(); ++i) {
            for (size_t j = i + 1; j < nb.size(); ++j) {
                const int64_t gap = nb[j].z_um - nb[i].z_um;
                if (gap >= horizon_um) break;  // at the target is not under it
                if (nb[i].seg_uid == nb[j].seg_uid) continue;  // one tube's own two ends

                std::array<int64_t, 4> key{
                    int64_t(std::min(nb[i].seg_uid, nb[j].seg_uid)),
                    int64_t(std::max(nb[i].seg_uid, nb[j].seg_uid)),
                    std::min(nb[i].z_um, nb[j].z_um),
                    std::max(nb[i].z_um, nb[j].z_um)};
                if (!seen_pairs.insert(key).second) continue;

                const bool  is_forced = shares_cell(nb[i].edge_idx, nb[j].edge_idx);
                Bucket     &b         = is_forced ? forced : avoidable;
                ++b.pairs;
                b.shortfall_um += horizon_um - gap;
                if (gap < tight_um) ++b.within_tight;
                if (gap == 0)       ++b.coincident;
                b.min_gap_um = std::min(b.min_gap_um, gap);
                if (!is_forced) ++cell_offences[cell];
            }
        }
    }

    BOOST_LOG_TRIVIAL(info) << label << " STAGGER: target "
        << std::fixed << std::setprecision(2) << horizon_um / 1000.0 << "mm | AVOIDABLE "
        << avoidable.pairs << " pairs under target ("
        << avoidable.within_tight << " under half, " << avoidable.coincident
        << " coincident), total shortfall "
        << avoidable.shortfall_um / 1000.0 << "mm, closest "
        << (avoidable.pairs ? avoidable.min_gap_um / 1000.0 : 0.0) << "mm"
        << " | FORCED (same-cell seams, unavoidable at full coverage) "
        << forced.pairs << " pairs, " << forced.coincident << " coincident";

    // Weak-plane census: cells with a boundary at each Z, over all boundaries, not split into
    // avoidable and forced.
    std::unordered_map<int64_t, int> plane_cells;
    std::vector<int64_t> zs;
    for (const auto &bounds_entry : cell_bounds) {
        const std::vector<Boundary> &bounds = bounds_entry.second;
        zs.clear();
        zs.reserve(bounds.size());
        for (const Boundary &b : bounds) zs.push_back(b.z_um);
        std::sort(zs.begin(), zs.end());
        zs.erase(std::unique(zs.begin(), zs.end()), zs.end());
        for (int64_t z : zs) ++plane_cells[z];
    }
    std::vector<std::pair<int, int64_t>> planes;
    planes.reserve(plane_cells.size());
    for (const auto &pc : plane_cells) planes.push_back({pc.second, pc.first});
    std::sort(planes.begin(), planes.end(), std::greater<>());

    const double total_cells = double(std::max<size_t>(1, m_cells.size()));
    if (!planes.empty())
        BOOST_LOG_TRIVIAL(info) << label << " PLANES: " << planes.size()
            << " distinct boundary Z, worst is Z=" << std::fixed << std::setprecision(3)
            << planes[0].second / 1000.0 << "mm shared by " << planes[0].first
            << "/" << m_cells.size() << " cells (" << std::setprecision(1)
            << 100.0 * double(planes[0].first) / total_cells << "%)";

    for (size_t i = 0; i < planes.size() && i < 10; ++i)
        BOOST_LOG_TRIVIAL(debug) << label << " PLANE #" << i << ": Z="
            << std::fixed << std::setprecision(3) << planes[i].second / 1000.0
            << "mm, " << planes[i].first << " cells ("
            << std::setprecision(1)
            << 100.0 * double(planes[i].first) / total_cells << "%)";

    std::vector<std::pair<size_t, CellId>> offenders;
    offenders.reserve(cell_offences.size());
    for (const auto &co : cell_offences) offenders.push_back({co.second, co.first});
    std::sort(offenders.begin(), offenders.end(),
              [](const auto &x, const auto &y) { return x.first > y.first; });
    for (size_t i = 0; i < offenders.size() && i < 10; ++i)
        BOOST_LOG_TRIVIAL(debug) << label << " CROWDED #" << i << ": cell ("
            << offenders[i].second.a << "," << offenders[i].second.b << ","
            << offenders[i].second.c << ") " << offenders[i].first
            << " avoidable pairs within " << std::fixed << std::setprecision(2)
            << horizon_um / 1000.0 << "mm";
}

} // namespace magma
} // namespace Slic3r
