/**
 * nisan_prg.h - Nisan's Pseudorandom Generator for Space-Bounded Computation
 * ==========================================================================
 * Declares data structures and functions for Nisan's PRG construction (1992).
 *
 * Knowledge: L3 (PairwiseHashFamily, SpacePRG), L5 (Nisan PRG algorithm),
 * L7 (BPL derandomization application).
 *
 * Refs: Nisan (1992), Arora-Barak Section 20.3.
 */

#ifndef NISAN_PRG_H
#define NISAN_PRG_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include "space_derand.h"

/* L3: Pairwise independent hash family (concrete instance).
 * Wraps PairwiseHashFamily with Nisan-specific parameters. */
typedef struct {
    size_t n;              /* Input bits */
    size_t m;              /* Output bits */
    size_t num_hashes;     /* Number of hash functions in this instance */
    PairwiseHashFamily *family; /* Underlying family */
} NisanHashDesc;

/* L5: Create a Nisan PRG for given space and time bounds.
 * Seed length: O(space * log time). Output length: time. */
SpacePRG *nisan_prg_create(size_t space, size_t output_len);

/* L5: Evaluate Nisan's PRG on a seed to produce pseudorandom output. */
bool nisan_prg_eval(const SpacePRG *prg, const bool *seed, size_t slen,
                     bool *out, size_t olen);

/* L5: Seed length required for Nisan's PRG. */
size_t nisan_prg_seed_len(size_t space, size_t output_len);

/* L5: Output length achievable from given seed length. */
size_t nisan_prg_output_len(size_t seed_len, size_t space);

/* L5: Verify statistical uniformity of PRG output. */
bool nisan_prg_verify_uniformity(const SpacePRG *prg, size_t trials);

/* L5: Free PRG resources. */
void nisan_prg_free(SpacePRG *prg);

/* L3: Create a pairwise independent hash family. */
PairwiseHashFamily *nisan_hash_create(size_t n, size_t m);

/* L3: Evaluate a hash function at index idx on input x. */
size_t nisan_hash_eval(const PairwiseHashFamily *h, size_t idx, size_t x);

/* L3: Verify pairwise independence property. */
bool nisan_hash_pairwise_independent(const PairwiseHashFamily *h, size_t trials);

/* L3: Free hash family. */
void nisan_hash_free(PairwiseHashFamily *h);

/* L7: Derandomize a BPL computation using PRG enumeration.
 * Implements BPL subseteq DSPACE(log^2 n) via Nisan's PRG. */
bool nisan_derandomize_bpl(size_t space, size_t time,
                            bool (*test_fn)(const bool*, size_t, void*),
                            void *ctx, bool *result);

/* L5: Error bound for given seed length in Nisan's construction. */
double nisan_error_for_seed_len(size_t space, size_t output_len,
                                 size_t seed_len);

/* L5: Optimal recursion depth for Nisan's tree construction. */
size_t nisan_optimal_depth(size_t output_len);

/* L5: Space used during PRG evaluation (O(space + log output_len)). */
size_t nisan_eval_space(size_t space, size_t output_len);

/* L8: Stretch factor of Nisan's PRG: output_len / seed_len.
 * For space S = O(log n): stretch = n^{Omega(1/log n)}. */
double nisan_stretch_factor(size_t space, size_t output_len);

/* L7: Application: derandomized algorithm for undirected s-t connectivity.
 * Uses Nisan's PRG instead of true random bits, maintaining correctness
 * with high probability. */
bool nisan_ustconn_derandomized(size_t n, const size_t *edges, size_t m,
                                 size_t s, size_t t, bool *connected);

/* L7: Application: derandomized bipartite matching verification.
 * Checks if a given matching is maximum using PRG-based random walks. */
bool nisan_matching_verify(size_t n_left, size_t n_right,
                            const size_t *edges, size_t m,
                            const size_t *matching, size_t match_size,
                            bool *is_maximum);

/* L8: Comparison: Nisan PRG vs truly random bits for space-bounded tests.
 * Statistical distance between PRG output distribution and uniform. */
double nisan_statistical_distance(size_t space, size_t output_len,
                                   size_t seed_len);

#endif /* NISAN_PRG_H */
