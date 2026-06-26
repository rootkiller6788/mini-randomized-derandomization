#ifndef DISP_CORE_H
#define DISP_CORE_H

/**
 * mini-randomness-extractors — Dispersers
 *
 * A disperser is a weaker variant of an extractor: instead of producing
 * a nearly-uniform output, it guarantees that the support of the output
 * is large (hits a large fraction of the output space).
 *
 * L1 Definitions: Disperser, oblivious disperser
 * L2 Core Concepts: hitting property vs extraction property
 * L4 Fundamental Laws: equivalence of dispersers to hitting set generators
 *
 * References:
 *   Sipser (1988) "Expanders, Randomness, or Time vs Space"
 *   Zuckerman (1996) "Simulating BPP Using a General Weak Random Source"
 *   Goldreich & Wigderson (1997) "Tiny Families of Functions..."
 */

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>

/*──────────────────────────────────────────────────────────────────────
 * L1: Disperser type
 *──────────────────────────────────────────────────────────────────────*/

/** A (k,ε)-disperser Disp : {0,1}^n × {0,1}^d → {0,1}^m.
 *
 *  For every distribution X over {0,1}^n with H∞(X) ≥ k,
 *  |Supp(Disp(X, U_d))| ≥ (1-ε) · 2^m.
 *
 *  That is: the output hits almost the entire range (large support),
 *  but does not need to be statistically close to uniform.
 *
 *  n = source bits, d = seed bits, m = output bits
 */
typedef struct {
    size_t n;        /**< Source length (bits) */
    size_t d;        /**< Seed length (bits) */
    size_t m;        /**< Output length (bits) */
} Disperser;

/** An oblivious disperser: a deterministic function O : {0,1}^n → {0,1}^m
 *  such that for every subset S ⊆ {0,1}^n of size ≥ 2^k,
 *  |O(S)| ≥ (1-ε) · 2^m.
 *
 *  This is a zero-seed disperser — no randomness needed beyond the source.
 */
typedef struct {
    size_t n;        /**< Source length (bits) */
    size_t d;        /**< Seed length (0 for truly oblivious) */
    size_t m;        /**< Output length (bits) */
    double eps;      /**< Error ε (fraction of range not hit) */
} ObliviousDisp;

/*──────────────────────────────────────────────────────────────────────
 * Core operations
 *──────────────────────────────────────────────────────────────────────*/

/** Create a disperser with given parameters.
 *  Validates: d ≤ n, m ≤ n, n ≤ 32 (for explicit construction).
 */
Disperser disp_create(size_t n, size_t d, size_t m);

/** Evaluate the disperser on a concrete source and seed.
 *  src: n bools representing the source bits
 *  src_len: must equal n
 *  seed: d bools representing the seed
 *  out: output buffer of m bools
 *  Returns true iff evaluation succeeds.
 */
bool disp_eval(const Disperser *d, const bool *src, size_t src_len,
               const bool *seed, bool *out);

/** Verify the hitting property of the disperser by random testing.
 *  Draws random subsets of size ≥ 2^k, checks that output support
 *  covers ≥ (1-ε) · 2^m of the range.
 *  k is derived from the construction parameters.
 */
bool disp_hits_all(const Disperser *d, size_t trials);

/*──────────────────────────────────────────────────────────────────────
 * Oblivious disperser
 *──────────────────────────────────────────────────────────────────────*/

/** Construct a zero-error oblivious disperser using algebraic
 *  (polynomial-based) construction.
 *
 *  Construction: View input x ∈ GF(2^n) as field element.
 *  Output: O(x) = (p_1(x), ..., p_m(x)) where p_i are degree-O(n/m)
 *  polynomials over GF(2^n).
 *
 *  This guarantees: for any K of size ≥ k, |O(K)| = 2^m.
 */
ObliviousDisp odisp_zero_error(size_t n, size_t d, size_t m);

/** Verify the oblivious disperser guarantee empirically.
 *  Enumerates all source positions (feasible for n ≤ 24) and
 *  checks the hitting property.
 */
bool odisp_verify(const ObliviousDisp *od, size_t trials);

/*──────────────────────────────────────────────────────────────────────
 * L4: Parameter bounds
 *──────────────────────────────────────────────────────────────────────*/

/** Minimum degree (2^d) required for given parameters.
 *  Based on the probabilistic method bound:
 *    D ≥ (n-k+log(1/ε)+O(1)) / (m - log(1/ε)) · O(1)
 */
size_t disp_degree(size_t n, double min_entropy, double error);

/** Maximum output bits achievable from given source.
 *  m ≤ k - 1 - log(1/ε)  (disperser bound, slightly better than extractor).
 */
size_t disp_output_bits(size_t n, double min_entropy);

/** Error probability bound for a random disperser.
 *  ε ≤ 2^{k} · 2^{d} · (2^m choose 2^m · (1-ε')) / 2^{n+d}
 *  Simplified: ε ≤ 2^{k + d - n - d + m} = 2^{k + m - n}.
 */
double disp_error_prob(size_t n, size_t d, size_t m, double k);

/** Construct a disperser from an extractor.
 *  Since every (k,ε)-extractor is also a (k,ε)-disperser,
 *  this is a trivial wrapper that demonstrates the relationship.
 *
 *  Theorem: If Ext is a (k,ε)-extractor, then Ext is a (k,ε)-disperser.
 *  Proof: Δ(Ext(X,U_d), U_m) ≤ ε  ⇒  |Supp(Ext(X,U_d))| ≥ (1-ε)·2^m.
 */
Disperser disp_from_extractor(size_t n, size_t d, size_t m);

/** Construct a disperser via bipartite expander graph.
 *  Uses the explicit expander construction from extractors
 *  (the "extractor graph" in extractor_constructions.h).
 *
 *  This is a non-trivial reduction: expander graphs with
 *  appropriate parameters yield dispersers.
 */
Disperser disp_from_expander(size_t n, size_t d, size_t m);

#endif /* DISP_CORE_H */
