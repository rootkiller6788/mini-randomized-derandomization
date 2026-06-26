#ifndef EXP_CORE_H
#define EXP_CORE_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>

/* ------------------------------------------------------------------
 * Core Data Types
 * ------------------------------------------------------------------ */

/** ExpanderGraph: d-regular expander on n vertices.
 *  adj: flat array of size n*d, adj[u*d + k] = k-th neighbor of u.
 *  lambda: cached second eigenvalue (0 if not computed). */
typedef struct { size_t n, d; size_t *adj; double lambda; } ExpanderGraph;

/** GraphVector: vector of length n over graph vertices. */
typedef struct { double *vec; size_t n; } GraphVector;

/* ------------------------------------------------------------------
 * Graph Lifecycle
 * ------------------------------------------------------------------ */

ExpanderGraph *exp_create(size_t n, size_t d);
void exp_free(ExpanderGraph *g);

/* ------------------------------------------------------------------
 * Edge Operations
 * ------------------------------------------------------------------ */

bool exp_add_edge(ExpanderGraph *g, size_t u, size_t v);
bool exp_has_edge(const ExpanderGraph *g, size_t u, size_t v);
size_t exp_neighbors(const ExpanderGraph *g, size_t u, size_t *neighbors);
size_t exp_degree(const ExpanderGraph *g, size_t u);

/* ------------------------------------------------------------------
 * Regularity & Connectivity
 * ------------------------------------------------------------------ */

bool exp_is_regular(const ExpanderGraph *g);
bool exp_is_connected(const ExpanderGraph *g);
bool exp_is_bipartite(const ExpanderGraph *g);

/* ------------------------------------------------------------------
 * Spectral Quantities
 * ------------------------------------------------------------------ */

double exp_second_eigenvalue(const ExpanderGraph *g);
double exp_spectral_gap(const ExpanderGraph *g);
double exp_rayleigh_quotient(const ExpanderGraph *g, const GraphVector *v);
double exp_collision_probability(const ExpanderGraph *g, size_t steps);
double exp_entropy_rate(const ExpanderGraph *g);

/* ------------------------------------------------------------------
 * Expander Verification
 * ------------------------------------------------------------------ */

bool exp_is_expander(const ExpanderGraph *g, double gamma);
bool exp_is_lambda_expander(const ExpanderGraph *g, double lambda);

/* ------------------------------------------------------------------
 * Cheeger Constant (Edge Expansion)
 * ------------------------------------------------------------------ */

double exp_cheeger_constant(const ExpanderGraph *g);
double exp_cheeger_lower_bound(const ExpanderGraph *g);
void exp_cheeger_bounds(const ExpanderGraph *g, double *lower, double *upper);
double exp_edge_expansion_exact(const ExpanderGraph *g);

/* ------------------------------------------------------------------
 * Vertex Expansion
 * ------------------------------------------------------------------ */

bool exp_verify_expansion(const ExpanderGraph *g, size_t S, double phi);
double exp_vertex_expansion_ratio(const ExpanderGraph *g, const char *subset);

/* ------------------------------------------------------------------
 * Random Walks & Mixing
 * ------------------------------------------------------------------ */

size_t exp_random_walk_mix_time(const ExpanderGraph *g, double eps);
size_t exp_mixing_time_precise(const ExpanderGraph *g, double eps);
size_t exp_random_walk_step(const ExpanderGraph *g, size_t current);
void exp_random_walk_distribution(const ExpanderGraph *g, size_t start,
                                   size_t steps, size_t samples, double *dist);
double exp_total_variation_distance(const double *p, size_t n);

/* ------------------------------------------------------------------
 * Graph Statistics
 * ------------------------------------------------------------------ */

size_t exp_num_edges(const ExpanderGraph *g);
size_t exp_diameter(const ExpanderGraph *g);
size_t exp_girth_lower_bound(const ExpanderGraph *g);

/* ------------------------------------------------------------------
 * GraphVector Operations
 * ------------------------------------------------------------------ */

GraphVector *gv_create(size_t n);
void gv_free(GraphVector *v);
double gv_norm(const GraphVector *v);
void gv_normalize(GraphVector *v);
double gv_dot(const GraphVector *a, const GraphVector *b);
void exp_apply_adjacency(const ExpanderGraph *g, const GraphVector *v,
                          GraphVector *result);

/* ------------------------------------------------------------------
 * Analysis & Display
 * ------------------------------------------------------------------ */

void exp_print_summary(const ExpanderGraph *g, FILE *fp);
double exp_alon_boppana_bound_static(size_t d);
double exp_compare_to_optimal(const ExpanderGraph *g);

/* ------------------------------------------------------------------
 * Spectral Theory (expander_spectral.c)
 * ------------------------------------------------------------------ */

double **exp_laplacian_matrix(const ExpanderGraph *g);
double **exp_normalized_laplacian(const ExpanderGraph *g);
void exp_laplacian_free(double **L, size_t n);
double exp_algebraic_connectivity(const ExpanderGraph *g);
double *exp_fiedler_vector(const ExpanderGraph *g, int iters);

/* Expander Mixing Lemma */
double exp_mixing_lemma_bound(const ExpanderGraph *g,
                               size_t S_size, size_t T_size);
size_t exp_count_edges_between(const ExpanderGraph *g,
                                 const char *S, const char *T);
double exp_verify_mixing_lemma(const ExpanderGraph *g, size_t max_size);

/* Independence & Chromatic Number Bounds */
double exp_independence_number_bound(const ExpanderGraph *g);
size_t exp_independence_number_exact(const ExpanderGraph *g);
double exp_chromatic_number_bound(const ExpanderGraph *g);

/* Full Spectrum (Jacobi method) */
int exp_full_spectrum(const ExpanderGraph *g, double *eigenvalues);

/* Trace Operations */
double exp_trace_adjacency(const ExpanderGraph *g);
double exp_trace_A2(const ExpanderGraph *g);
double exp_trace_A3(const ExpanderGraph *g);
size_t exp_num_triangles(const ExpanderGraph *g);

/* Heat Kernel */
double exp_heat_kernel_diagonal(const ExpanderGraph *g, double t);
double exp_heat_kernel_mixing_time(const ExpanderGraph *g, double eps);

/* Spectral Partitioning */
bool exp_spectral_partition(const ExpanderGraph *g, char *part, int iters);
size_t exp_partition_cut_size(const ExpanderGraph *g, const char *part);
double exp_partition_balance(const ExpanderGraph *g, const char *part);

#endif /* EXP_CORE_H */
