/**
 * logspace_derand.h - Logspace Derandomization: USTCONN and Branching Programs
 * ==========================================================================
 * Declares logspace-specific data structures and functions for
 * derandomizing space-bounded computation.
 *
 * Knowledge: L1 (LogspaceFn, BranchingProgram), L4 (Reingold Theorem),
 * L5 (logspace derandomization), L6 (USTCONN canonical problem).
 *
 * Refs: Reingold (2008), Arora-Barak Ch.4, 20.
 */

#ifndef LOGSPACE_DERAND_H
#define LOGSPACE_DERAND_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include "space_derand.h"

/* L1: Logspace-computable function representation.
 * A function f: {0,1}^n -> {0,1} computable in O(log n) space. */
typedef struct {
    size_t n;                      /* Input length */
    bool (*fn)(const bool*, size_t);/* Function pointer */
} LogspaceFn;

/* L5: Logspace derandomization via PRG enumeration. */
bool sd_logspace_derandomize(const BranchingProgram *bp, size_t n,
                              bool *result);

/* L4: USTCONN verification (undirected graph check). */
bool sd_verify_undirected(size_t n, const size_t *edges, size_t m);

/* L5: Component counting for undirected graphs. */
size_t sd_count_components_undirected(size_t n, const size_t *edges, size_t m);

/* L5: Build expander from a base graph via zig-zag iterations. */
RegularGraph *sd_build_expander_from_graph(const RegularGraph *G, size_t rounds);

/* L4: Reingold seed length for logspace derandomization. */
size_t sd_reingold_seed_len(size_t n, size_t d, size_t layers);

/* L5: Diameter bound computation for branching programs. */
bool sd_bp_diameter_bound_compute(const BranchingProgram *bp, size_t *diameter);

/* L5: Mixing time estimate for random walks on BPs. */
size_t sd_bp_mixing_time(const BranchingProgram *bp);

/* L5: Check if a branching program is deterministic. */
bool sd_bp_is_deterministic(const BranchingProgram *bp);

/* L6: USTCONN - undirected s-t connectivity in logspace (Reingold 2008). */
bool sd_ustconn_logspace(size_t n, const size_t *edges, size_t m,
                          size_t s, size_t t, bool *connected);

/* L6: STCONN - directed s-t connectivity in NL (via Savitch). */
bool sd_stconn_nl(size_t n, const size_t *edges, size_t m,
                   size_t s, size_t t, bool *connected);

/* L6: PATH - check if a path exists in a directed graph of bounded depth. */
bool sd_path_exists_bounded(const ConfigGraph *cg, size_t s, size_t t,
                             size_t max_depth, bool *found);

/* L3: Create a logspace-computable rotation map for a d-regular graph. */
bool sd_rotation_map_create(size_t n, size_t d, const size_t *edges,
                             size_t *rotation);

/* L5: Simulate a Turing machine step in logspace.
 * Given current config and transition table, compute next config. */
bool sd_tm_step(const SpaceTM *tm, const char *input, size_t ilen,
                 const TMConfig *cur, TMConfig *next);

/* L7: Application: derandomized graph connectivity for network routing.
 * Uses Reingold's algorithm for reliable connectivity testing in
 * space-constrained environments (embedded systems, sensor networks). */
bool sd_network_connectivity_check(size_t n_nodes, const size_t *links,
                                    size_t num_links, bool *all_connected);

/* L7: Application: model checking safety properties in logspace.
 * Verifies that no unsafe state is reachable from initial state. */
bool sd_safety_property_check(const ConfigGraph *sys, const size_t *unsafe,
                               size_t num_unsafe, bool *safe);

/* L8: Derandomized squaring of graphs for improved expansion.
 * Advanced technique for boosting spectral gap without randomness. */
RegularGraph *sd_derandomized_square(const RegularGraph *G, size_t steps);

#endif /* LOGSPACE_DERAND_H */
