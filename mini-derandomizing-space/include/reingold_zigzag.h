/**
 * reingold_zigzag.h - Reingold Zig-Zag Product and Expander Constructions
 * ==========================================================================
 * Declares graph product operations and expander analysis functions for
 * Reingold's logspace USTCONN algorithm (2008).
 *
 * Knowledge: L3 (RegularGraph, SpectralData), L5 (zig-zag, power, tensor,
 * replacement products), L8 (expander mixing lemma, spectral analysis).
 *
 * Refs: Reingold-Vadhan-Wigderson (2002), Hoory-Linial-Wigderson (2006).
 */

#ifndef REINGOLD_ZIGZAG_H
#define REINGOLD_ZIGZAG_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include "space_derand.h"

/* L3: Create a d-regular graph from an edge list (rotation map). */
RegularGraph *rz_create_from_edges(size_t n, size_t d, const size_t *edges);

/* L3: Create an n-cycle (2-regular, simplest expander). */
RegularGraph *rz_create_cycle(size_t n);

/* L3: Create a complete graph K_n ((n-1)-regular, perfect expander). */
RegularGraph *rz_create_complete(size_t n);

/* L5: Zig-zag product G (R) H.  |V| = |V_G|*|V_H|, degree = d_H^2.
 * Core operation in Reingold's construction.
 * Ref: Reingold-Vadhan-Wigderson (2002) Theorem 3.2. */
RegularGraph *rz_zigzag_product(const RegularGraph *G, const RegularGraph *H);

/* L5: Graph power G^t. Connect vertices within distance t.
 * Degree = d^t. Used to boost spectral gap. */
RegularGraph *rz_power_graph(const RegularGraph *G, size_t t);

/* L5: Tensor product G x H. Eigenvalues multiply.
 * Ref: Hoory-Linial-Wigderson Proposition 4.3. */
RegularGraph *rz_tensor_product(const RegularGraph *G, const RegularGraph *H);

/* L5: Replacement product G (r) H. Degree = d_H + 1.
 * Ref: Reingold-Vadhan-Wigderson Section 3.3. */
RegularGraph *rz_replacement_product(const RegularGraph *G, const RegularGraph *H);

/* L5: Check if G is an (n,d,lambda)-expander. */
bool rz_is_expander(const RegularGraph *G, double lambda);

/* L5: Compute spectral gap gamma = 1 - lambda_2. */
double rz_spectral_gap(const RegularGraph *G);

/* L5: Free a regular graph. */
void rz_free(RegularGraph *G);

/* L4: Build logspace connectivity network (Reingold 2008 Section 3).
 * Iteratively applies squaring + zig-zag to produce an expander. */
RegularGraph *rz_logspace_connectivity_network(size_t n);

/* L4: Verify graph connectivity via DFS. */
bool rz_verify_connectivity(const RegularGraph *G);

/* L5: Shortest path in unweighted regular graph (BFS). */
size_t rz_shortest_path(const RegularGraph *G, size_t s, size_t t);

/* L4: Verify Reingold's theorem for given parameters. */
bool rz_reingold_theorem_check(size_t n, size_t d);

/* L3: Compute normalized adjacency matrix (row-major). */
bool rz_adjacency_matrix(const RegularGraph *G, double *mat);

/* L5: Compute all eigenvalues via QR/Power iteration. */
bool rz_compute_all_eigenvalues(const RegularGraph *G, double *evals,
                                  size_t max_n);

/* L8: Expander mixing lemma verification.
 * |E(S,T) - d|S||T|/n| <= lambda * d * sqrt(|S||T|).
 * Ref: Hoory-Linial-Wigderson Lemma 2.4. */
bool rz_mixing_lemma_check(const RegularGraph *G,
                            const size_t *S, size_t s_size,
                            const size_t *T, size_t t_size,
                            double *edge_count, double *bound);

/* L8: Alon-Boppana bound: for any infinite family of d-regular graphs,
 * liminf lambda_2 >= 2*sqrt(d-1)/d. This verifies the bound. */
double rz_alon_boppana_bound(size_t d);

/* L8: Ramanujan graphs: graphs achieving lambda <= 2*sqrt(d-1)/d.
 * Check if a given graph is Ramanujan. */
bool rz_is_ramanujan(const RegularGraph *G);

/* L5: Compute the diameter of a regular graph. */
size_t rz_diameter(const RegularGraph *G);

/* L5: Random walk on regular graph for t steps.
 * Returns final vertex distribution (probability per vertex). */
bool rz_random_walk_distribution(const RegularGraph *G, size_t start,
                                  size_t steps, double *distribution);

/* L7: Application: derandomized routing in expander-based networks.
 * Uses zig-zag product to construct low-degree, high-expansion
 * networks for reliable packet routing. */
bool rz_routing_network_build(size_t n_nodes, size_t target_degree,
                               RegularGraph **network);

/* L7: Application: error-correcting codes from expander graphs.
 * Expander-based LDPC codes achieving near-Shannon capacity. */
bool rz_expander_code_check(const RegularGraph *G, size_t *min_distance);

/* L9: Research frontier: high-dimensional expanders.
 * Generalize spectral expansion to simplicial complexes.
 * Placeholder for future exploration. */
bool rz_high_dimensional_expander_test(const RegularGraph *G, size_t dim);

#endif /* REINGOLD_ZIGZAG_H */
