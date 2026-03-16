#include "MagmaGreedyWarmStart.hpp"

#include <boost/log/trivial.hpp>

#include <algorithm>
#include <chrono>
#include <queue>
#include <set>
#include <vector>

namespace Slic3r {
namespace magma {

// ============================================================================
// CellConsumed — tracks consumed Z-ranges per cell (sorted, non-overlapping)
// ============================================================================

namespace {

struct CellConsumed {
    // Sorted non-overlapping intervals (start_um, end_um)
    std::vector<std::pair<int64_t, int64_t>> intervals;

    bool overlaps(int64_t start_um, int64_t end_um) const {
        // Binary search for first interval that could overlap
        auto it = std::lower_bound(intervals.begin(), intervals.end(),
            std::make_pair(start_um, start_um),
            [](const std::pair<int64_t, int64_t> &a,
               const std::pair<int64_t, int64_t> &b) {
                return a.second <= b.first;
            });
        return it != intervals.end() && it->first < end_um;
    }

    void add(int64_t start_um, int64_t end_um) {
        // Find merge range: all existing intervals that overlap or touch [start, end]
        auto lo = std::lower_bound(intervals.begin(), intervals.end(),
            std::make_pair(start_um, start_um),
            [](const std::pair<int64_t, int64_t> &a,
               const std::pair<int64_t, int64_t> &b) {
                return a.second < b.first;
            });
        auto hi = lo;
        while (hi != intervals.end() && hi->first <= end_um) {
            start_um = std::min(start_um, hi->first);
            end_um   = std::max(end_um, hi->second);
            ++hi;
        }
        auto pos = intervals.erase(lo, hi);
        intervals.insert(pos, {start_um, end_um});
    }
};

// ============================================================================
// CellLayerScore — priority queue entry (min-heap on score)
// ============================================================================

struct CellLayerScore {
    TriangleCell cell;
    int          layer;
    int64_t      layer_bottom_um;
    int64_t      layer_top_um;
    double       score; // lower = more constrained = higher priority

    bool operator>(const CellLayerScore &o) const { return score > o.score; }
};

// Find the run in an edge that contains a given layer. Returns nullptr if none.
const Run *find_run_containing(const EdgeData &ed, int layer, int64_t bottom_um, int64_t top_um)
{
    for (const Run &r : ed.runs) {
        if (layer >= r.start_layer && layer <= r.end_layer &&
            bottom_um >= r.start_um && top_um <= r.end_um)
            return &r;
    }
    return nullptr;
}

// Count unconsumed layers in a run for a given cell
int count_unconsumed_layers(const Run &run,
                            const CellConsumed &consumed_a,
                            const CellConsumed &consumed_b,
                            const MicronTables &um)
{
    int count = 0;
    for (int L = run.start_layer; L <= run.end_layer; ++L) {
        int64_t bot = um.bottom_um[L];
        int64_t top = um.top_um[L];
        if (!consumed_a.overlaps(bot, top) && !consumed_b.overlaps(bot, top))
            ++count;
    }
    return count;
}

} // anonymous namespace

// ============================================================================
// greedy_warm_start
// ============================================================================

void greedy_warm_start(
    const std::unordered_map<TriangleCell, CellPresence, TriangleCellHash> &cells,
    const std::vector<EdgeData>                                             &edges,
    const std::unordered_map<TriangleCell, std::vector<size_t>, TriangleCellHash> &cell_edges,
    const MicronTables                                                      &um,
    int64_t min_h_um,
    int64_t max_h_um,
    std::vector<std::vector<CommittedSegment>>                              &committed)
{
    auto t_start = std::chrono::high_resolution_clock::now();

    // Consumed intervals per cell
    std::unordered_map<TriangleCell, CellConsumed, TriangleCellHash> consumed;

    using MinHeap = std::priority_queue<CellLayerScore, std::vector<CellLayerScore>,
                                         std::greater<CellLayerScore>>;

    // Periodic re-scoring interval. As tubes are assigned, the heap ordering
    // becomes stale (cell×layers that were "easy" may now be constrained).
    // Rebuilding the heap periodically corrects the priority ordering so the
    // most-constrained-first heuristic stays effective. Scaled by model size:
    // num_edges roughly tracks complexity, /3 gives ~3 re-scores during the
    // main assignment wave. Floor of 200 prevents churn on tiny models.
    const int rescore_every = std::max(200, static_cast<int>(edges.size()) / 3);

    // ------------------------------------------------------------------
    // Build heap: score unconsumed cell×layers
    // ------------------------------------------------------------------
    // Score = sum of achievable tube heights across all unconsumed neighbors.
    // Lower score = fewer/shorter options = more constrained = higher priority.

    auto build_heap = [&]() -> MinHeap {
        MinHeap heap;
        for (const auto &[cell, presence] : cells) {
            auto ce_it = cell_edges.find(cell);
            if (ce_it == cell_edges.end()) continue;

            for (int L = presence.first_layer; L <= presence.last_layer; ++L) {
                if (!presence.present(L)) continue;

                int64_t layer_bot = um.bottom_um[L];
                int64_t layer_top = um.top_um[L];

                if (consumed[cell].overlaps(layer_bot, layer_top)) continue;

                double score = 0.0;
                bool fillable = false;

                for (size_t ei : ce_it->second) {
                    const EdgeData &ed = edges[ei];
                    const Run *run = find_run_containing(ed, L, layer_bot, layer_top);
                    if (!run) continue;

                    const TriangleCell &neighbor =
                        (ed.edge.a == cell) ? ed.edge.b : ed.edge.a;
                    if (consumed[neighbor].overlaps(layer_bot, layer_top)) continue;

                    int64_t tube_start = layer_bot;
                    int64_t tube_end   = layer_top;
                    const CellConsumed &cons_c = consumed[cell];
                    const CellConsumed &cons_n = consumed[neighbor];

                    for (int lo = L - 1; lo >= run->start_layer; --lo) {
                        int64_t bot = um.bottom_um[lo];
                        if (tube_end - bot > max_h_um) break;
                        if (cons_c.overlaps(bot, um.top_um[lo]) ||
                            cons_n.overlaps(bot, um.top_um[lo])) break;
                        tube_start = bot;
                    }
                    for (int hi = L + 1; hi <= run->end_layer; ++hi) {
                        int64_t top = um.top_um[hi];
                        if (top - tube_start > max_h_um) break;
                        if (cons_c.overlaps(um.bottom_um[hi], top) ||
                            cons_n.overlaps(um.bottom_um[hi], top)) break;
                        tube_end = top;
                    }

                    int64_t potential = tube_end - tube_start;
                    if (potential >= min_h_um) {
                        fillable = true;
                        score += double(potential);
                    }
                }

                if (fillable)
                    heap.push({cell, L, layer_bot, layer_top, score});
            }
        }
        return heap;
    };

    // ------------------------------------------------------------------
    // Greedy assignment with periodic re-scoring
    // ------------------------------------------------------------------

    int total_assigned = 0, total_skipped = 0, rescores = 0;
    int since_rescore = rescore_every; // trigger initial build

    MinHeap heap;

    for (;;) {
        // Rebuild heap when stale (or on first iteration)
        if (since_rescore >= rescore_every) {
            heap = build_heap();
            ++rescores;
            since_rescore = 0;
            if (heap.empty()) break;
        }

        if (heap.empty()) break;

        CellLayerScore entry = heap.top();
        heap.pop();

        // Skip if consumed since last scoring
        if (consumed[entry.cell].overlaps(entry.layer_bottom_um, entry.layer_top_um)) {
            ++total_skipped;
            continue;
        }

        // Find the most constrained unconsumed neighbor at this Z
        auto ce_it = cell_edges.find(entry.cell);
        if (ce_it == cell_edges.end()) continue;

        size_t best_edge_idx = SIZE_MAX;
        const Run *best_run  = nullptr;
        int best_neighbor_free = INT_MAX;
        TriangleCell best_neighbor;

        for (size_t ei : ce_it->second) {
            const EdgeData &ed = edges[ei];
            const TriangleCell &neighbor =
                (ed.edge.a == entry.cell) ? ed.edge.b : ed.edge.a;

            auto pres_it = cells.find(neighbor);
            if (pres_it == cells.end() || !pres_it->second.present(entry.layer))
                continue;
            if (consumed[neighbor].overlaps(entry.layer_bottom_um, entry.layer_top_um))
                continue;

            const Run *run = find_run_containing(ed, entry.layer,
                                                  entry.layer_bottom_um, entry.layer_top_um);
            if (!run) continue;

            int free_layers = count_unconsumed_layers(*run,
                consumed[entry.cell], consumed[neighbor], um);

            if (free_layers < best_neighbor_free) {
                best_neighbor_free = free_layers;
                best_edge_idx      = ei;
                best_run           = run;
                best_neighbor      = neighbor;
            }
        }

        if (best_edge_idx == SIZE_MAX) continue;

        // Expand the longest valid tube containing this layer
        int64_t tube_start = entry.layer_bottom_um;
        int64_t tube_end   = entry.layer_top_um;
        const CellConsumed &cons_a = consumed[entry.cell];
        const CellConsumed &cons_b = consumed[best_neighbor];

        for (int lo = entry.layer - 1; lo >= best_run->start_layer; --lo) {
            int64_t bot = um.bottom_um[lo];
            int64_t top = um.top_um[lo];
            if (tube_end - bot > max_h_um) break;
            if (cons_a.overlaps(bot, top) || cons_b.overlaps(bot, top)) break;
            tube_start = bot;
        }

        for (int hi = entry.layer + 1; hi <= best_run->end_layer; ++hi) {
            int64_t bot = um.bottom_um[hi];
            int64_t top = um.top_um[hi];
            if (top - tube_start > max_h_um) break;
            if (cons_a.overlaps(bot, top) || cons_b.overlaps(bot, top)) break;
            tube_end = top;
        }

        if (tube_end - tube_start < min_h_um) continue;

        committed[best_edge_idx].push_back({tube_start, tube_end});
        consumed[entry.cell].add(tube_start, tube_end);
        consumed[best_neighbor].add(tube_start, tube_end);
        ++total_assigned;
        ++since_rescore;
    }

    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - t_start).count();
    BOOST_LOG_TRIVIAL(info) << "MagmaGreedy: " << total_assigned << " tubes assigned, "
        << total_skipped << " skipped, " << rescores << " rescores, " << ms << "ms";
}

} // namespace magma
} // namespace Slic3r
