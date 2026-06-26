#ifndef EXP_APPS_H
#define EXP_APPS_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include "expander_core.h"

/* ------------------------------------------------------------------
 * Data Types
 * ------------------------------------------------------------------ */

/** Random walk configuration: graph + number of steps. */
typedef struct { const ExpanderGraph *g; size_t steps; } RandomWalk;

/* ------------------------------------------------------------------
 * Error Reduction & Derandomization
 * ------------------------------------------------------------------ */

/** Amplify RP/BPP error probability using expander random walk. */
bool exp_error_reduction(const ExpanderGraph *g, size_t trials,
                          double base_err, double *amplified_err);

/** Amplify an RP algorithm (one-sided error) via expander walk. */
bool exp_amplify_rp_algorithm(const ExpanderGraph *g,
                               bool (*rp)(void), size_t trials,
                               double *error);

/** Derandomized sampling: generate approximately uniform samples
 *  using expander random walk (saves random bits). */
bool exp_derandomized_sampling(const ExpanderGraph *g, size_t samples,
                                size_t *indices);

/** Quality of derandomized sampling: TV distance to uniform. */
double exp_sampling_quality(const size_t *samples, size_t num_samples,
                             size_t n);

/** Random bits saved by expander walk vs independent sampling. */
double exp_random_bits_saved(size_t k, size_t n, size_t d);

/* ------------------------------------------------------------------
 * SL = L Connectivity
 * ------------------------------------------------------------------ */

/** Test s-t connectivity using spectral expander properties. */
bool exp_sl_connectivity(const ExpanderGraph *g, size_t s, size_t t,
                          bool *connected);

/** Estimate shortest path length between s and t (BFS or spectral). */
size_t exp_sl_path_length(const ExpanderGraph *g, size_t s, size_t t);

/* ------------------------------------------------------------------
 * Sorting Networks
 * ------------------------------------------------------------------ */

/** Generate comparator pairs for sorting network from expander topology. */
bool exp_sorting_network_from_expander(size_t n, size_t *comparators_out,
                                        size_t *num_comparators);

/* ------------------------------------------------------------------
 * Chernoff Bounds on Expanders
 * ------------------------------------------------------------------ */

/** Chernoff-style concentration bound for expander random walk
 *  (Gillman 1998). */
double exp_chernoff_bound_on_expander(const ExpanderGraph *g, size_t n,
                                       double mu, double delta);

/** Standard Chernoff bound for independent samples (comparison). */
double exp_chernoff_bound_independent(size_t n, double mu, double delta);

/* ------------------------------------------------------------------
 * Expander Codes
 * ------------------------------------------------------------------ */

/** Encode binary message using expander-based error-correcting code. */
bool exp_code_encode(const char *message, size_t m,
                      const ExpanderGraph *g, char *codeword_out);

/** Iterative bit-flip decoding for expander codes. */
bool exp_code_flip_decode(char *codeword, size_t m,
                           const ExpanderGraph *g, size_t max_iters);

/* ------------------------------------------------------------------
 * Pseudorandom Generators & Extractors
 * ------------------------------------------------------------------ */

/** Generate pseudorandom sequence via expander random walk. */
bool exp_prg_generate(const ExpanderGraph *g, size_t start,
                       size_t length, size_t *output);

/** Simulate an extractor using expander walk (source + seed -> output). */
bool exp_extractor_graph(const ExpanderGraph *g, const char *source,
                          size_t seed, size_t length, char *output);

/* ------------------------------------------------------------------
 * Lossless Expanders & Hitting Times
 * ------------------------------------------------------------------ */

/** Check if graph is a lossless expander (heuristic via sampling). */
bool exp_is_lossless_expander(const ExpanderGraph *g, size_t k,
                               double eps, size_t num_samples);

/** Expected hitting time from s to t via random walk simulation. */
double exp_random_walk_hitting_time(const ExpanderGraph *g, size_t s,
                                     size_t t, size_t num_simulations,
                                     size_t max_steps);

/* ------------------------------------------------------------------
 * Load Balancing
 * ------------------------------------------------------------------ */

/** Distribute loads across processors using expander diffusion. */
void exp_load_balance(const double *load, const ExpanderGraph *g,
                       double *new_load, size_t rounds);

/** Maximum load imbalance ratio (0 = perfect balance). */
double exp_load_imbalance(const double *load, size_t n);

#endif /* EXP_APPS_H */
