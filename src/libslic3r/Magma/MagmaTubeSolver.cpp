#include "MagmaTubeSolver.hpp"
#include "MagmaGreedyWarmStart.hpp"

#include <boost/log/trivial.hpp>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>

// OR-Tools CP-SAT solver
#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/sat/cp_model_checker.h"

namespace Slic3r {
namespace magma {

// ============================================================================
// Constructor
// ============================================================================

MagmaTubeSolver::MagmaTubeSolver(
    const std::unordered_map<TriangleCell, CellPresence, TriangleCellHash> &cells,
    const std::vector<LayerData> &layer_data,
    double min_tube_height_mm,
    double max_tube_height_mm,
    int    num_layers,
    double stagger_period_mm,
    MagmaTubeSolverMode mode,
    double solver_timeout_sec,
    int    stagger_tolerance_pct)
    : m_cells(cells)
    , m_layer_data(layer_data)
    , m_min_h_mm(min_tube_height_mm)
    , m_max_h_mm(max_tube_height_mm)
    , m_num_layers(num_layers)
    , m_stagger_period_mm(stagger_period_mm)
    , m_mode(mode)
    , m_timeout_sec(solver_timeout_sec)
    , m_stagger_tolerance_pct(stagger_tolerance_pct)
{}

// ============================================================================
// solve — 3-pass driver
// ============================================================================

void MagmaTubeSolver::solve(
    std::vector<UTubePair> &out_pairs,
    std::unordered_map<TriangleCell, std::vector<int>, TriangleCellHash> &out_cell_pair_index,
    ProgressFn progress_fn,
    ThrowIfCanceled throw_if_canceled)
{
    auto t_start = std::chrono::high_resolution_clock::now();

    build_micron_tables();
    build_edges();
    if (throw_if_canceled) throw_if_canceled();

    // Greedy warm start: fast heuristic that populates m_committed with a good
    // initial solution. CP-SAT then refines it via warm start hints.
    {
        const int64_t min_h_um = llround(m_min_h_mm * 1000.0);
        const int64_t max_h_um = llround(m_max_h_mm * 1000.0);
        greedy_warm_start(m_cells, m_edges, m_cell_edges, m_um,
                          min_h_um, max_h_um, m_committed,
                          &m_cell_difficulty);
        validate_committed(m_edges, m_edge_index, m_committed, m_cells,
                           m_um, m_layer_data, m_min_h_mm, m_max_h_mm,
                           m_num_layers, "GREEDY");
    }
    if (throw_if_canceled) throw_if_canceled();

    // Z window and stride, derived from max tube height:
    //   window  = 4 × max_h_layers  — room for 3-4 stacked tubes
    //   overlap = 2 × max_h_layers  — every tube fully in ≥2 Z levels
    //   stride  = 2 × max_h_layers  — window minus overlap
    // Uses smallest layer height for conservative layer count.
    int z_stride;
    {
        double min_lh = m_layer_data.empty() ? 0.2 : m_layer_data[0].height;
        for (const auto &ld : m_layer_data)
            if (ld.height > 0 && ld.height < min_lh)
                min_lh = ld.height;
        int max_h_layers = std::max(1, static_cast<int>(std::ceil(m_max_h_mm / min_lh)));
        // Full-Z: window covers entire object height.
        m_z_window = m_num_layers;
        z_stride   = std::max(1, m_num_layers);
    }

    // Z extent of actual edge data (not all layers may have cells)
    int max_edge_layer = 0;
    for (const auto &[cell, presence] : m_cells)
        max_edge_layer = std::max(max_edge_layer, presence.last_layer);

    int num_z_levels = 0;
    for (int z = 0; z <= max_edge_layer; z += z_stride)
        ++num_z_levels;

    BOOST_LOG_TRIVIAL(info) << "MagmaTubeSolver: " << m_edges.size() << " edges"
        << ", Z_window=" << m_z_window << " layers"
        << ", Z_stride=" << z_stride
        << ", R=" << R
        << ", " << num_z_levels << " Z levels"
        << ", max_layer=" << max_edge_layer;

    // CP-SAT refinement (Refined mode only)
    if (m_mode == MagmaTubeSolverMode::Refined) {
        int current_z = 0;
        if (progress_fn) progress_fn(0, num_z_levels);
        for (int z_off = 0; z_off <= max_edge_layer; z_off += z_stride) {
            solve_pass(0,     0,     z_off);
            if (throw_if_canceled) throw_if_canceled();
            solve_pass(R / 2, R / 2, z_off);
            if (throw_if_canceled) throw_if_canceled();
            if (progress_fn) progress_fn(++current_z, num_z_levels);
        }

        validate_committed(m_edges, m_edge_index, m_committed, m_cells,
                           m_um, m_layer_data, m_min_h_mm, m_max_h_mm,
                           m_num_layers, "CPSAT");
    }

    extract_results(out_pairs, out_cell_pair_index);

    auto t_end = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();
    BOOST_LOG_TRIVIAL(info) << "MagmaTubeSolver: "
        << (m_mode == MagmaTubeSolverMode::Refined ? "refined" : "basic")
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

    for (int L = 0; L < n; ++L)
        m_um.top_um[L] = llround(m_layer_data[L].print_z * 1000.0);

    // Layer 0 bottom from actual bottom_z; subsequent layers contiguous by construction
    m_um.bottom_um[0] = llround(m_layer_data[0].bottom_z() * 1000.0);
    for (int L = 1; L < n; ++L)
        m_um.bottom_um[L] = m_um.top_um[L - 1];

    // Reverse maps (same integers — guaranteed to match)
    for (int L = 0; L < n; ++L) {
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
        for (const TriangleCell &nbr : cell.neighbors()) {
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

            if (!ed.runs.empty()) {
                size_t idx = m_edges.size();
                m_edge_index[edge] = idx;
                m_edges.push_back(std::move(ed));
            }
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

void MagmaTubeSolver::build_blocks(int off_a, int off_b, int off_z,
                                    std::vector<Block> &out) const
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

    // cell -> XY block key
    auto cell_block = [&](const TriangleCell &c) -> XYKey {
        // Use floor division with offset
        int shifted_a = c.a + off_a;
        int shifted_b = c.b + off_b;
        int bx = (shifted_a >= 0) ? shifted_a / R : (shifted_a - R + 1) / R;
        int by = (shifted_b >= 0) ? shifted_b / R : (shifted_b - R + 1) / R;
        return {bx, by};
    };

    // Group cells by XY block
    std::unordered_map<XYKey, std::vector<TriangleCell>, XYKeyHash> xy_groups;
    for (const auto &[cell, _] : m_cells)
        xy_groups[cell_block(cell)].push_back(cell);

    // Single Z slice: off_z to off_z + z_window - 1
    // (Z levels are iterated by the caller in solve())
    int z_start = std::max(0, off_z);
    int z_end   = std::min(off_z + m_z_window - 1, m_num_layers - 1);
    if (z_start >= m_num_layers || z_start > z_end) return;

    int64_t z_start_um = m_um.bottom_um[z_start];
    int64_t z_end_um   = m_um.top_um[z_end];

    // Build blocks: one per XY group at this Z slice
    for (const auto &[xy_key, cells] : xy_groups) {
        std::unordered_set<TriangleCell, TriangleCellHash> cell_set(cells.begin(), cells.end());

        Block block;
        block.cells = cell_set;
        block.z_start_layer = z_start;
        block.z_end_layer   = z_end;
        block.z_start_um    = z_start_um;
        block.z_end_um      = z_end_um;

        // Collect edges via cell reverse lookup (avoids scanning all edges)
        std::unordered_set<size_t> seen_edges;
        for (const TriangleCell &cell : cells) {
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

        if (!block.edge_indices.empty())
            out.push_back(std::move(block));
    }

    BOOST_LOG_TRIVIAL(debug) << "MagmaTubeSolver: pass offset(" << off_a << ","
        << off_b << ") -> " << out.size() << " blocks";
}

// ============================================================================
// solve_block — build and solve CP-SAT model for one block
// ============================================================================

BlockResult MagmaTubeSolver::solve_block(const Block &block) const
{
    using namespace operations_research::sat;

    const int64_t min_h_um = llround(m_min_h_mm * 1000.0);
    const int64_t max_h_um = llround(m_max_h_mm * 1000.0);

    // ------------------------------------------------------------------
    // Phase-offset stagger grid
    // ------------------------------------------------------------------
    // WHY 3 GRIDS: triangle 3-coloring guarantees adjacent edges (which share
    // a cell) always use different phase grids. Since each edge's SharedEdge
    // type (Horizontal/Col60/Diag120) maps deterministically to a phase,
    // stagger is structural — it requires zero explicit pairwise constraints.
    //
    // WHY DOMAIN RESTRICTION (not objective penalties): the original design
    // used W_GRID_DIST and W_DODGE soft objectives, but domain restriction is
    // superior for two reasons:
    //   1. Smaller domains → tighter LP relaxation → faster CP-SAT solving.
    //      Boundaries go from ~10-30 positions (all layers) to ~3-8 (grid +
    //      endpoints), and the O(n²) feasible-size computation shrinks too.
    //   2. Dodge penalties spread injections across many Z levels (penalizing
    //      nearby pairs pushes boundaries apart everywhere). Domain restriction
    //      concentrates injections on a small set of grid-aligned Z levels,
    //      which is preferable for practical injection logistics.
    //
    // FILL-RATE GUARANTEE: run endpoints are always included in the domain
    // (see lines below), so the maximum-coverage solution is always feasible
    // regardless of grid spacing. Grid points provide additional stagger
    // options within runs, but never remove the ability to achieve full fill.
    const int64_t period_um = llround(m_stagger_period_mm * 1000.0);
    const int64_t z0_um = m_um.bottom_um.empty() ? 0 : m_um.bottom_um[0];

    // Pre-compute the 3 phase grid snap sets (sorted allowed boundary positions)
    std::array<std::unordered_set<int64_t>, 3> phase_snap;
    if (period_um > 0) {
        // Collect all layer boundaries
        std::set<int64_t> all_bounds;
        for (int L = 0; L < m_num_layers; ++L) {
            all_bounds.insert(m_um.bottom_um[L]);
            all_bounds.insert(m_um.top_um[L]);
        }
        int64_t z_max = m_um.top_um[m_num_layers - 1];
        int64_t phase_offset = period_um / 3;

        for (int phase = 0; phase < 3; ++phase) {
            int64_t offset = phase * phase_offset;
            for (int64_t g = z0_um + offset; g <= z_max + period_um; g += period_um) {
                // Snap to nearest actual layer boundary
                auto it = all_bounds.lower_bound(g);
                int64_t best = -1, best_dist = INT64_MAX;
                if (it != all_bounds.end() && (*it - g) < best_dist)
                    { best_dist = *it - g; best = *it; }
                if (it != all_bounds.begin())
                    { --it; if ((g - *it) < best_dist) best = *it; }
                if (best >= 0) phase_snap[phase].insert(best);
            }
        }

        BOOST_LOG_TRIVIAL(debug) << "MagmaTubeSolver: phase grids: "
            << phase_snap[0].size() << "/" << phase_snap[1].size() << "/"
            << phase_snap[2].size() << " positions (period=" << period_um/1000.0 << "mm)";
    }

    CpModelBuilder model;

    // Track intervals per cell for NoOverlap
    std::unordered_map<TriangleCell, std::vector<IntervalVar>, TriangleCellHash> cell_intervals;

    // Track segment variables for objective and result extraction
    struct SegVars {
        BoolVar     active;
        IntVar      start;
        IntVar      end;
        IntVar      size;
        IntervalVar interval;
        size_t      edge_idx;
        int64_t     run_start_um;  // which run this slot belongs to
        int64_t     run_end_um;
        operations_research::Domain contrib_dom; // {0} ∪ feasible sizes
        std::vector<int64_t> boundaries;              // boundary values for this run
        std::unordered_set<int64_t> on_grid;          // grid + run endpoint positions
    };
    std::vector<SegVars> all_segments;

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

            // K = enough slots to fill this run with average-sized tubes, plus
            // 1 for stagger splitting. Also floor at the greedy tube count for
            // this edge (so the solver can always represent the greedy solution).
            // With full-Z, runs can span the entire object (a "run" is an
            // uninterrupted cell-pair presence range, broken by model edges or
            // mid-Z constrictions).
            // K = greedy tube count in this run + 1 (room for one stagger split).
            // Only count greedy tubes within this run's effective Z range.
            int64_t run_start_um = m_um.bottom_um[eff_start];
            int64_t run_end_um   = m_um.top_um[eff_end];
            int K_from_greedy = 0;
            for (const CommittedSegment &seg : m_committed[ei])
                if (seg.start_um >= run_start_um && seg.end_um <= run_end_um)
                    ++K_from_greedy;
            int K = std::max(1, K_from_greedy * 2 + 1);

            // Unified layer boundary list for this run.
            // When stagger grid is active, filter to the edge's phase grid.
            // Phase is determined by SharedEdge type (3-coloring).
            std::vector<int64_t> all_boundaries;
            all_boundaries.push_back(m_um.bottom_um[eff_start]);
            for (int L = eff_start; L <= eff_end; ++L)
                all_boundaries.push_back(m_um.top_um[L]);

            // Domain = {run endpoints} ∪ {grid-aligned layer boundaries}.
            //
            // Run endpoints are ALWAYS included — this is the fill-rate guarantee.
            // Without them, a run whose endpoints don't fall on the grid would
            // have no valid start/end position, losing coverage entirely.
            //
            // Grid points are additional options within runs. The solver picks
            // them when they enable better cross-edge packing via NoOverlap
            // (stagger "falls out" of the optimization). The resulting domain is
            // typically 3-8 values instead of 10-30 (all_boundaries), which
            // reduces the O(n²) feasible-size computation below and tightens
            // CP-SAT's LP relaxation for faster solving.
            // ----------------------------------------------------------
            // Difficulty-based domain selection.
            //
            // For each run, compute the average difficulty of both cells
            // across the run's layer range. Easy runs (low difficulty)
            // get grid-restricted domains for stagger. Hard runs (high
            // difficulty) get full domains to preserve coverage.
            //
            // Difficulty comes from greedy's initial unconstrained scoring:
            //   0 = easiest (3 neighbors at max_h)
            //   3 × max_h_um = hardest (no fillable neighbors)
            //
            // Threshold: 30% of max_possible = restrict when ≥70% of
            // theoretical max flexibility is available.
            // ----------------------------------------------------------
            const int64_t max_possible = 3 * max_h_um;
            const int64_t difficulty_threshold = max_possible * 3 / 10; // 30%

            int64_t run_difficulty_sum = 0;
            int run_difficulty_count = 0;
            auto lookup_diff = [&](const TriangleCell &cell, int layer) -> int64_t {
                auto it = m_cell_difficulty.find(cell);
                if (it == m_cell_difficulty.end()) return max_possible;
                auto pres_it = m_cells.find(cell);
                if (pres_it == m_cells.end()) return max_possible;
                int idx = layer - pres_it->second.first_layer;
                if (idx < 0 || idx >= static_cast<int>(it->second.size())) return max_possible;
                return it->second[idx];
            };
            for (int L = eff_start; L <= eff_end; ++L) {
                run_difficulty_sum += lookup_diff(ed.edge.a, L) + lookup_diff(ed.edge.b, L);
                run_difficulty_count += 2;
            }
            int64_t avg_difficulty = (run_difficulty_count > 0)
                ? run_difficulty_sum / run_difficulty_count : max_possible;
            if (period_um <= 0) {
                BOOST_LOG_TRIVIAL(error) << "MagmaTubeSolver: stagger period must be > 0";
                return BlockResult{};
            }

            bool restrict_domain = (avg_difficulty < difficulty_threshold);

            std::vector<int64_t> boundaries;
            std::unordered_set<int64_t> on_grid; // positions not penalized by W_OFFGRID
            if (restrict_domain) {
                // Easy run: grid + run endpoints + greedy hints.
                // Stagger-friendly domain with off-grid penalty for greedy positions.
                int phase = static_cast<int>(shared_edge(ed.edge.a, ed.edge.b));
                const auto &snap = phase_snap[phase];
                std::set<int64_t> bset;
                bset.insert(all_boundaries.front());
                bset.insert(all_boundaries.back());
                on_grid.insert(all_boundaries.front());
                on_grid.insert(all_boundaries.back());
                for (int64_t b : all_boundaries)
                    if (snap.count(b)) {
                        bset.insert(b);
                        on_grid.insert(b);
                    }
                // Greedy hint positions: valid domain but penalized (off-grid)
                int64_t run_lo = all_boundaries.front();
                int64_t run_hi = all_boundaries.back();
                for (const CommittedSegment &seg : m_committed[ei]) {
                    if (seg.start_um >= run_lo && seg.start_um <= run_hi)
                        bset.insert(seg.start_um);
                    if (seg.end_um >= run_lo && seg.end_um <= run_hi)
                        bset.insert(seg.end_um);
                }
                boundaries.assign(bset.begin(), bset.end());
            } else {
                // Hard run: full domain. All layer boundaries available,
                // all treated as on-grid (no off-grid penalty).
                boundaries = std::move(all_boundaries);
                for (int64_t b : boundaries)
                    on_grid.insert(b);
            }
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
                    // When both active, k must come after k-1 (gap of at least 1 um)
                    model.AddLessThan(all_segments.back().end, start)
                        .OnlyEnforceIf({prev_active, active});
                }

                int64_t rs = boundaries.front();
                int64_t re = boundaries.back();
                all_segments.push_back({active, start, end, size, interval,
                                        ei, rs, re, contrib_dom, boundaries, on_grid});
                prev_active = active;
                has_prev = true;
            }
        }
    }

    // ------------------------------------------------------------------
    // 2. Add frozen intervals from committed tubes outside the block
    // ------------------------------------------------------------------
    // Two kinds of frozen intervals prevent the solver from creating tubes
    // that overlap with committed segments it can't modify:
    //
    // XY boundary frozen: edges that cross the XY block boundary have one
    // cell inside and one outside. The committed tubes on these edges are
    // controlled by whichever block contains both cells. We freeze them on
    // the inside cell to prevent NoOverlap violations.
    //
    // KEY INVARIANT: a block only erases/replaces segments it can fully see
    // (both cells in XY). With full-Z, all decision-edge segments are fully
    // inside the Z range, so only XY boundary freezing is needed.

    // (A) XY boundary: scan edges crossing the block boundary
    std::unordered_set<size_t> seen_boundary_edges;
    for (const TriangleCell &cell : block.cells) {
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

            const TriangleCell &inside_cell = a_in ? ed.edge.a : ed.edge.b;

            for (const CommittedSegment &seg : m_committed[ei]) {
                // Clip to block Z range
                int64_t clip_start = std::max(seg.start_um, block.z_start_um);
                int64_t clip_end   = std::min(seg.end_um, block.z_end_um);
                if (clip_start >= clip_end)
                    continue;

                IntervalVar frozen = model.NewFixedSizeIntervalVar(
                    clip_start, clip_end - clip_start);
                cell_intervals[inside_cell].push_back(frozen);
            }
        }
    }

    // (B) Z boundary: decision-edge segments that extend outside the block's
    // Z range are frozen. Segments fully inside Z range become warm start
    // hints — the solver re-optimizes them.
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
        }
    }

    // Early out: no decision variables means nothing to solve
    if (all_segments.empty())
        return BlockResult{};

    // ------------------------------------------------------------------
    // 3. NoOverlap per cell
    // ------------------------------------------------------------------
    for (auto &[cell, intervals] : cell_intervals) {
        if (intervals.size() > 1)
            model.AddNoOverlap(intervals);
    }

    // ------------------------------------------------------------------
    // 4. Objective: coverage + excess tube budget + off-grid penalty
    //
    //   Tier 1 — Coverage:  W_COVERAGE × Σ{contrib}       (dominant)
    //   Tier 2 — Excess:    W_EXCESS × max(0, active - budget)
    //   Tier 3 — Off-grid:  W_OFFGRID × Σ{is_offgrid}     (tiebreaker)
    //
    // Stagger comes from domain restriction (phase-offset grids).
    // The excess budget allows the solver to add tubes for stagger up
    // to greedy_count × (1 + tolerance/100) without penalty. Beyond
    // that, each extra tube costs W_EXCESS — preventing fragmentation
    // while allowing stagger splits.
    //
    // Off-grid penalty nudges boundaries toward grid positions when
    // coverage is equal. Greedy hint positions are in the domain so
    // hints are valid, but the solver prefers grid positions.
    // ------------------------------------------------------------------

    constexpr int64_t W_COVERAGE = 1000000;
    constexpr int64_t W_EXCESS  = 100;
    constexpr int64_t W_OFFGRID = 1;

    LinearExpr objective;

    for (const auto &seg : all_segments) {
        IntVar contrib = model.NewIntVar(seg.contrib_dom);
        model.AddEquality(contrib, seg.size).OnlyEnforceIf(seg.active);
        model.AddEquality(contrib, 0).OnlyEnforceIf(seg.active.Not());
        objective += W_COVERAGE * contrib;
    }

    // Excess tube budget: free tubes up to greedy count × (1 + tolerance%).
    // Only tubes beyond the budget are penalized.
    {
        int greedy_count = 0;
        for (size_t ei : block.edge_indices)
            greedy_count += static_cast<int>(m_committed[ei].size());

        int budget = greedy_count * (100 + m_stagger_tolerance_pct) / 100;
        if (budget > 0) {
            LinearExpr total_active;
            for (const auto &seg : all_segments)
                total_active += seg.active;

            IntVar excess = model.NewIntVar(
                operations_research::Domain(0, static_cast<int>(all_segments.size())));
            model.AddGreaterOrEqual(excess, total_active - budget);
            objective -= W_EXCESS * excess;
        }
    }

    // Off-grid penalty: nudge boundaries toward grid positions.
    // Each boundary at a greedy (non-grid) position costs 1.
    // Max total per block: ~160 (80 segs × 2 boundaries). Negligible vs coverage.
    if (period_um > 0) {
        for (const auto &seg : all_segments) {
            auto make_offgrid = [&](IntVar boundary_var,
                                    const std::vector<int64_t> &bounds,
                                    const std::unordered_set<int64_t> &on_grid) -> IntVar {
                bool has_offgrid = false;
                for (int64_t b : bounds)
                    if (!on_grid.count(b)) { has_offgrid = true; break; }
                if (!has_offgrid)
                    return model.NewConstant(0);

                IntVar is_off = model.NewIntVar(operations_research::Domain(0, 1));
                auto table = model.AddAllowedAssignments({boundary_var, is_off});
                for (int64_t b : bounds)
                    table.AddTuple({b, on_grid.count(b) ? 0 : 1});
                return is_off;
            };

            IntVar start_off = make_offgrid(seg.start, seg.boundaries, seg.on_grid);
            IntVar end_off   = make_offgrid(seg.end, seg.boundaries, seg.on_grid);
            objective -= W_OFFGRID * start_off;
            objective -= W_OFFGRID * end_off;
        }
    }

    model.Maximize(objective);

    // ------------------------------------------------------------------
    // 5. Complete solution hint (warm start)
    // ------------------------------------------------------------------
    // Greedy segments become CP-SAT hints, letting the solver start from
    // ~80% coverage instead of searching from scratch. This is critical for
    // converging within the per-block timeout — without hints, blocks often
    // time out at UNKNOWN with poor or no solutions.
    //
    // Hinting inactive slots (below) gives a COMPLETE initial solution. With
    // all variables assigned, CP-SAT can immediately verify feasibility and
    // begin improvement search, rather than spending time constructing an
    // initial feasible solution that may be worse than greedy.
    {
        std::vector<bool> hinted(all_segments.size(), false);
        int hint_matched = 0, hint_missed = 0, total_greedy = 0;

        // For each edge in this block, try to match each greedy segment to
        // a solver slot in the same run. Greedy segments are sorted by start.
        for (size_t ei : block.edge_indices) {
            // Collect greedy segments for this edge within block Z range
            std::vector<const CommittedSegment *> greedy_segs;
            for (const CommittedSegment &seg : m_committed[ei])
                if (seg.start_um >= block.z_start_um && seg.end_um <= block.z_end_um)
                    greedy_segs.push_back(&seg);
            std::sort(greedy_segs.begin(), greedy_segs.end(),
                      [](const CommittedSegment *a, const CommittedSegment *b) {
                          return a->start_um < b->start_um;
                      });
            total_greedy += static_cast<int>(greedy_segs.size());

            // Match each greedy segment to the first available slot in its run
            for (const CommittedSegment *cs : greedy_segs) {
                bool matched = false;
                for (size_t si = 0; si < all_segments.size(); ++si) {
                    if (hinted[si]) continue;
                    const auto &sv = all_segments[si];
                    if (sv.edge_idx != ei) continue;
                    if (cs->start_um >= sv.run_start_um && cs->end_um <= sv.run_end_um) {
                        model.AddHint(sv.active, true);
                        model.AddHint(sv.start, cs->start_um);
                        model.AddHint(sv.end, cs->end_um);
                        hinted[si] = true;
                        matched = true;
                        ++hint_matched;
                        break;
                    }
                }
                if (!matched)
                    ++hint_missed;
            }
        }

        // Hint remaining slots as inactive (complete initial solution)
        for (size_t si = 0; si < all_segments.size(); ++si)
            if (!hinted[si])
                model.AddHint(all_segments[si].active, false);

        if (hint_missed > 0)
            BOOST_LOG_TRIVIAL(warning) << "MagmaTubeSolver: " << hint_missed
                << "/" << total_greedy << " greedy hints unmatched ("
                << hint_matched << " matched, "
                << (all_segments.size() - hint_matched) << " inactive)";
    }

    // ------------------------------------------------------------------
    // 6. Solve
    // ------------------------------------------------------------------
    SatParameters params;
    params.set_max_time_in_seconds(m_timeout_sec);
    params.set_num_workers(CPSAT_WORKERS);
    // NOTE: linearization_level and relative_gap_limit removed — defaults
    // perform better with our discrete domains. gap_limit was causing
    // medium blocks to stop early, losing 1-2 tubes per block.
    // NOTE: repair_hint and hint_conflict_limit removed — triggered
    // CP-SAT fixed_search crash (integer_search.cc:1217)

    // Validate model before solving (catches structural bugs)
    auto proto = model.Build();
    std::string validation_error = ValidateCpModel(proto);
    if (!validation_error.empty()) {
        BOOST_LOG_TRIVIAL(error) << "MagmaTubeSolver: MODEL INVALID: "
            << validation_error.substr(0, 200);
    }

    operations_research::sat::Model sat_model;
    sat_model.Add(NewSatParameters(params));

    // Cancellation: checked between passes via check_canceled().
    // Mid-solve abort not wired — TimeLimit::RegisterExternalBooleanAsLimit
    // conflicts with NewSatParameters (creates TimeLimit without params).

    CpSolverResponse response = SolveCpModel(proto, &sat_model);

    // ------------------------------------------------------------------
    // 7. Extract results
    // ------------------------------------------------------------------
    BlockResult result;
    if (response.status() == CpSolverStatus::OPTIMAL ||
        response.status() == CpSolverStatus::FEASIBLE) {
        result.solved = true;
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
        int64_t greedy_cov = 0;
        for (size_t ei : block.edge_indices)
            for (const CommittedSegment &cs : m_committed[ei])
                greedy_cov += cs.end_um - cs.start_um;

        BOOST_LOG_TRIVIAL(debug) << "MagmaTubeSolver: block solved "
            << (response.status() == CpSolverStatus::OPTIMAL ? "OPTIMAL" : "FEASIBLE")
            << ", obj=" << response.objective_value()
            << ", gap=" << (response.best_objective_bound() > 0
                ? (1.0 - response.objective_value() / response.best_objective_bound()) * 100.0
                : 0.0) << "%"
            << ", segs=" << result.segments.size()
            << "/" << all_segments.size()
            << ", wall=" << std::fixed << std::setprecision(2)
            << response.wall_time() << "s"
            << " | cov=" << solver_cov/1000.0
            << "mm greedy=" << greedy_cov/1000.0 << "mm";
        if (solver_cov < greedy_cov)
            BOOST_LOG_TRIVIAL(warning) << "  COVERAGE DROP: lost "
                << (greedy_cov - solver_cov)/1000.0 << "mm"
                << " (solver " << result.segments.size()
                << " segs vs greedy on " << block.edge_indices.size() << " edges)";
    } else {
        // Count frozen intervals per cell for diagnostics
        int total_frozen = 0;
        int max_frozen_per_cell = 0;
        for (const auto &[cell, intervals] : cell_intervals) {
            // Count non-decision intervals (frozen = total - decision)
            int n = static_cast<int>(intervals.size());
            if (n > max_frozen_per_cell) max_frozen_per_cell = n;
            total_frozen += n;
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
            << ", total_intervals=" << total_frozen
            << ", max_per_cell=" << max_frozen_per_cell;
    }

    return result;
}

// ============================================================================
// solve_pass
// ============================================================================

void MagmaTubeSolver::solve_pass(int off_a, int off_b, int off_z)
{
    auto t_start = std::chrono::high_resolution_clock::now();

    std::vector<Block> blocks;
    build_blocks(off_a, off_b, off_z, blocks);

    // Solve blocks in parallel — each is independent (no shared cells)
    std::vector<BlockResult> results(blocks.size());

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, blocks.size()),
        [&](const tbb::blocked_range<size_t> &range) {
            for (size_t i = range.begin(); i < range.end(); ++i)
                results[i] = solve_block(blocks[i]);
        }
    );

    // Merge results sequentially (no race — parallel_for is done)
    commit_results(blocks, results);

    auto t_end = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

    int total_segs = 0;
    for (const auto &br : results)
        total_segs += static_cast<int>(br.segments.size());

    BOOST_LOG_TRIVIAL(debug) << "MagmaTubeSolver: pass(" << off_a << "," << off_b
        << "): " << blocks.size() << " blocks, "
        << total_segs << " segments, " << ms << "ms";
}

// ============================================================================
// commit_results
// ============================================================================

void MagmaTubeSolver::commit_results(const std::vector<Block> &blocks,
                                      const std::vector<BlockResult> &results)
{
    // Only update edges from blocks that solved successfully.
    // If a block failed (INFEASIBLE/UNKNOWN), preserve prior committed segments
    // rather than destroying them — don't make things worse.
    //
    // For successful blocks: only remove segments FULLY INSIDE the block's Z
    // range. Segments extending outside Z belong to other Z levels and must be
    // preserved. This is the key invariant: a block only modifies what it can
    // fully see.
    for (size_t bi = 0; bi < blocks.size(); ++bi) {
        if (!results[bi].solved) {
            ++m_unknown_blocks;
            continue;
        }

        const Block &block = blocks[bi];
        for (size_t ei : block.edge_indices) {
            auto &segs = m_committed[ei];
            segs.erase(std::remove_if(segs.begin(), segs.end(),
                [&](const CommittedSegment &s) {
                    return s.start_um >= block.z_start_um &&
                           s.end_um   <= block.z_end_um;
                }), segs.end());
        }
        for (const auto &[edge_idx, seg] : results[bi].segments)
            m_committed[edge_idx].push_back(seg);
    }
}

// ============================================================================
// extract_results — convert micron-space assignments to UTubePair format
// ============================================================================

void MagmaTubeSolver::extract_results(
    std::vector<UTubePair> &out_pairs,
    std::unordered_map<TriangleCell, std::vector<int>, TriangleCellHash> &out_index) const
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

    validate_committed(m_edges, m_edge_index, m_committed, m_cells,
                       m_um, m_layer_data, m_min_h_mm, m_max_h_mm,
                       m_num_layers, "FINAL");
}

// ============================================================================
// validate_committed — reusable validation for committed segments
// ============================================================================

ValidationResult validate_committed(
    const std::vector<EdgeData>                                            &edges,
    const std::unordered_map<CellEdge, size_t, CellEdgeHash>              &edge_index,
    const std::vector<std::vector<CommittedSegment>>                       &committed,
    const std::unordered_map<TriangleCell, CellPresence, TriangleCellHash> &cells,
    const MicronTables                                                     &um,
    const std::vector<LayerData>                                           &layer_data,
    double min_h_mm, double max_h_mm, int num_layers,
    const char *label)
{
    ValidationResult result;

    // Build per-cell segment list for overlap checking
    std::unordered_map<TriangleCell, std::vector<std::pair<int, int>>,
                       TriangleCellHash> cell_segments; // (start_layer, end_layer)

    int total_segs = 0;
    for (size_t ei = 0; ei < edges.size(); ++ei) {
        const EdgeData &ed = edges[ei];
        for (const CommittedSegment &seg : committed[ei]) {
            ++total_segs;
            auto start_it = um.bottom_to_layer.find(seg.start_um);
            auto end_it   = um.top_to_layer.find(seg.end_um);
            if (start_it == um.bottom_to_layer.end() ||
                end_it   == um.top_to_layer.end()) {
                if (++result.bad_range <= 3)
                    BOOST_LOG_TRIVIAL(warning) << label << ": invalid micron boundary "
                        << seg.start_um << "-" << seg.end_um;
                continue;
            }

            int sL = start_it->second;
            int eL = end_it->second;

            // 1. Valid layer range
            if (sL < 0 || eL >= num_layers || sL > eL) {
                if (++result.bad_range <= 3)
                    BOOST_LOG_TRIVIAL(warning) << label << ": tube invalid range L"
                        << sL << "-" << eL;
                continue;
            }

            // 2. Height bounds
            double h_mm = layer_data[eL].print_z
                        - (layer_data[sL].print_z - layer_data[sL].height);
            if (h_mm < min_h_mm - 0.01) {
                if (++result.bad_short <= 3)
                    BOOST_LOG_TRIVIAL(warning) << label << ": tube too short: "
                        << h_mm << "mm (min=" << min_h_mm << ") L" << sL << "-" << eL;
            }
            if (h_mm > max_h_mm + 0.01) {
                if (++result.bad_long <= 3)
                    BOOST_LOG_TRIVIAL(warning) << label << ": tube too long: "
                        << h_mm << "mm (max=" << max_h_mm << ") L" << sL << "-" << eL;
            }

            // 3. Edge exists
            CellEdge ce(ed.edge.a, ed.edge.b);
            if (edge_index.find(ce) == edge_index.end()) {
                if (++result.bad_edge <= 3)
                    BOOST_LOG_TRIVIAL(warning) << label << ": cells are not a valid edge";
            }

            // 4. Both cells present at every layer
            auto it_a = cells.find(ed.edge.a);
            auto it_b = cells.find(ed.edge.b);
            if (it_a != cells.end() && it_b != cells.end()) {
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

    // 5. Per-cell overlap
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

    // 6. Per-edge overlap
    for (size_t ei = 0; ei < edges.size(); ++ei) {
        const auto &segs = committed[ei];
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
            << result.bad_range << " bad-range, " << result.bad_edge << " bad-edge, "
            << result.bad_presence << " not-present, "
            << result.overlap_cell << " cell-overlap, "
            << result.overlap_edge << " edge-overlap)";
    else
        BOOST_LOG_TRIVIAL(info) << label << ": all " << total_segs << " segments OK";

    // Coverage summary
    {
        double total_presence_um = 0, total_covered_um = 0;
        int cells_below_50 = 0, cells_below_25 = 0, cells_zero = 0;
        for (const auto &[cell, presence] : cells) {
            int64_t presence_um = um.top_um[presence.last_layer]
                                - um.bottom_um[presence.first_layer];
            if (presence_um <= 0) continue;

            int64_t covered_um = 0;
            auto it = cell_segments.find(cell);
            if (it != cell_segments.end()) {
                for (const auto &[sL, eL] : it->second) {
                    covered_um += um.top_um[eL] - um.bottom_um[sL];
                }
                // Each segment counted twice (once per cell), but we're iterating
                // per-cell so each cell sees its own segments. However segments
                // appear in cell_segments for BOTH cells of the edge, so each
                // cell correctly sees all its tubes.
            }
            double pct = double(covered_um) / double(presence_um);
            total_presence_um += presence_um;
            total_covered_um  += covered_um;
            if (covered_um == 0) ++cells_zero;
            else if (pct < 0.25) ++cells_below_25;
            else if (pct < 0.50) ++cells_below_50;
        }
        double overall_pct = total_presence_um > 0
            ? total_covered_um / total_presence_um * 100.0 : 0.0;
        BOOST_LOG_TRIVIAL(info) << label << " COVERAGE: " << std::fixed << std::setprecision(1)
            << overall_pct << "% overall, "
            << cells_zero << " cells unfilled, "
            << cells_below_25 << " cells <25%, "
            << cells_below_50 << " cells <50%"
            << " (covered=" << int64_t(total_covered_um) << "um"
            << " presence=" << int64_t(total_presence_um) << "um)";
    }

    return result;
}

} // namespace magma
} // namespace Slic3r
