#ifndef EXP_CONSTR_H
#define EXP_CONSTR_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include "expander_core.h"

/* ------------------------------------------------------------------
 * Construction Parameter Types
 * ------------------------------------------------------------------ */

/** Parameters for LPS Ramanujan construction: primes p, q. */
typedef struct { size_t p, q; } RamanujanParams;

/** Parameters for zig-zag product: n vertices, d degree, rotation map. */
typedef struct { size_t n, d; size_t *rotations; } ZigZagParams;

/* ------------------------------------------------------------------
 * Explicit Constructions
 * ------------------------------------------------------------------ */

/** Margulis-Gabber-Galil 8-regular expander on m^2 vertices. */
ExpanderGraph *exp_margulis_gabber_galil(size_t n);

/** LPS Ramanujan graph: degree p+1, approximately q(q^2-1)/2 vertices. */
ExpanderGraph *exp_ramanujan_lps(size_t p, size_t q);

/** Random d-regular graph via configuration model. */
ExpanderGraph *exp_random_regular(size_t n, size_t d);

/** Cayley graph on Z_n with given generator set. */
ExpanderGraph *exp_cayley_graph(size_t n, const size_t *generators,
                                 size_t num_gen);

/* ------------------------------------------------------------------
 * Graph Products & Operations
 * ------------------------------------------------------------------ */

/** Zig-zag product G (z) H. Returns the result ExpanderGraph directly. */
ExpanderGraph *exp_zigzag_product(const ExpanderGraph *G,
                                   const ExpanderGraph *H);

/** Replacement product G (r) H. */
ExpanderGraph *exp_replacement_product(const ExpanderGraph *G,
                                        const ExpanderGraph *H);

/** Tensor product G1 x G2. */
ExpanderGraph *exp_tensor_product(const ExpanderGraph *G1,
                                   const ExpanderGraph *G2);

/** k-th graph power G^k. */
ExpanderGraph *exp_graph_power(const ExpanderGraph *G, size_t k);

/* ------------------------------------------------------------------
 * Ramanujan & Spectral Properties
 * ------------------------------------------------------------------ */

/** Check if graph satisfies the Ramanujan bound lambda_2 <= 2*sqrt(d-1). */
bool exp_is_ramanujan(const ExpanderGraph *g);

/** Static check: does given lambda satisfy Ramanujan bound? */
bool exp_is_ramanujan_bound(double lambda2, size_t d);

/** Alon-Boppana bound: liminf lambda_2 >= 2*sqrt(d-1). */
double exp_alon_boppana_bound(size_t d);

/** Nilli (finite-n) bound: lambda_2 >= 2*sqrt(d-1) * (1 - O(1/log_k(n))). */
double exp_nilli_bound(size_t n, size_t d);

/** Expected lambda_2 for large random d-regular graphs. */
double exp_expected_second_eigenvalue(size_t d);

/** Heuristic: spectral gap >= 10% of degree => good expander. */
bool exp_is_good_expander(const ExpanderGraph *g);

/** Search for best Cayley expander generators (small n). */
double exp_find_best_cayley_expander(size_t n, size_t d,
                                      size_t *best_generators);

#endif /* EXP_CONSTR_H */
