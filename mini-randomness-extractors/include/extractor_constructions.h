#ifndef EXT_CONSTR_H
#define EXT_CONSTR_H

/**
 * mini-randomness-extractors — Extractor Constructions
 *
 * L5 Algorithms: Leftover Hash Lemma, Trevisan, Raz, Ta-Shma,
 *                 Carter-Wegman hashing, two-source extractors
 * L8 Advanced: extractor graphs, design-based constructions
 *
 * References:
 *   Impagliazzo, Levin, Luby (1989) "Pseudo-random Generation from
 *     One-way Functions" — Leftover Hash Lemma
 *   Trevisan (2001) "Extractors and Pseudorandom Generators"
 *   Raz (2005) "Extractors with Weak Random Seeds"
 *   Ta-Shma (1996) "On Extracting Randomness From Weak Random Sources"
 *   Chor & Goldreich (1988) "Unbiased Bits from Sources of Weak
 *     Randomness and a Communication Complexity Puzzle"
 *   Bourgain (2005) "More on the Sum-Product Phenomenon..."
 *   Carter & Wegman (1979) "Universal Classes of Hash Functions"
 */

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>

#include "extractor_core.h"

/*──────────────────────────────────────────────────────────────────────
 * L3: Universal Hash Functions — Carter-Wegman construction
 *──────────────────────────────────────────────────────────────────────*/

/** A family of universal hash functions H: {0,1}^n → {0,1}^m.
 *  Carter-Wegman: h_{a,b}(x) = ((a·x + b) mod p) mod 2^m
 *  where p is a prime > 2^n, a ∈ [1,p-1], b ∈ [0,p-1].
 *
 *  For pairwise independent hashing (stronger):
 *    h_{A,b}(x) = (A·x + b) over GF(2^m)
 */
typedef struct {
    size_t n;        /**< Input length (bits) */
    size_t m;        /**< Output length (bits) */
    size_t p;        /**< Prime modulus p > 2^n */
    size_t *a;       /**< Multiplier array (could be matrix for PI) */
    size_t *b;       /**< Offset array */
    size_t num_a;    /**< Number of a parameters */
    size_t num_b;    /**< Number of b parameters */
} UniversalHash;

/** Construct a Carter-Wegman universal hash family.
 *  H = { h_{a,b}(x) = ((a·x + b) mod p) mod 2^m | a∈[1,p-1], b∈[0,p-1] }
 *  Number of functions: (p-1) * p.
 */
UniversalHash uh_create_carter_wegman(size_t n, size_t m);

/** Evaluate hash function h_{a,b}(x) given parameters.
 *  x: input value (interpreted as integer from n bits)
 *  a_idx, b_idx: parameter indices (0-based)
 *  Returns h(x) as size_t (m bits).
 */
size_t uh_eval(const UniversalHash *uh, size_t a_idx, size_t b_idx,
               size_t x);

/** Evaluate with a single seed index (encodes both a and b).
 *  Useful as the extractor function interface.
 */
size_t uh_eval_seeded(const UniversalHash *uh, size_t seed, size_t x);

/** Verify the universal property: for any distinct x≠y,
 *  Pr_h[h(x)=h(y)] ≤ 1/2^m.
 */
bool uh_check_universal(const UniversalHash *uh, size_t trials);

/** Get total number of hash functions in the family.
 */
size_t uh_num_functions(const UniversalHash *uh);

/** Free universal hash family resources.
 */
void uh_free(UniversalHash *uh);

/** Create a pairwise independent hash family over GF(2^m).
 *  h_{A,b}(x) = A·x + b  (matrix-vector multiply over GF(2)).
 *  This is 2-universal: Pr[h(x)=a ∧ h(y)=b] = 1/2^{2m} for x≠y.
 */
UniversalHash uh_create_pairwise_independent(size_t n, size_t m);

/*──────────────────────────────────────────────────────────────────────
 * L4+L5: Extractor Constructions
 *──────────────────────────────────────────────────────────────────────*/

/** Leftover Hash Lemma (Impagliazzo-Levin-Luby) extractor.
 *
 *  Theorem (LHL): Let H be a universal hash family from {0,1}^n to
 *  {0,1}^m. For any random variable X over {0,1}^n with H∞(X) ≥ k,
 *  Δ((H, H(X)), (H, U_m)) ≤ 2^{-(k-m)/2}.
 *
 *  Returns a SeededExt where:
 *    - n is source length, m is output length
 *    - seed selects a hash function from the family
 *    - ε = 2^{-(k-m)/2} for source min-entropy k
 */
SeededExt ext_leftover_hash(size_t n, size_t m);

/** Trevisan extractor (2001): extractor from Nisan-Wigderson PRG.
 *
 *  Uses: (1) a weak combinatorial design, (2) a "hard" predicate
 *  (here: parity function), (3) the NW reconstruction paradigm.
 *
 *  For source length n, output m, error ε:
 *    seed length d = O(log n + log(1/ε))
 *    Can extract m = k^Ω(1) bits from source min-entropy k.
 *
 *  This implementation uses a truncated Hadamard code as the predicate
 *  and a simple algebraic design.
 */
SeededExt ext_trevisan(size_t n, size_t m, double eps);

/** Raz extractor (2005): extractor with weak random seed.
 *
 *  Handles the case where even the seed has imperfect randomness.
 *  Composition: inner extractor (seed → "good" seed) → outer extractor.
 *
 *  This implementation uses the block-and-condense paradigm.
 */
SeededExt ext_raz(size_t n, size_t m, double eps);

/** Ta-Shma extractor (1996): extractor based on error-correcting codes.
 *
 *  Key idea: encode the source using a binary error-correcting code,
 *  then output a random subset of the codeword bits.
 *
 *  Uses a Reed-Solomon-like concatenated code in this implementation.
 *  Achieves near-optimal entropy loss: m ≈ k - O(log(1/ε)).
 */
SeededExt ext_ta_shma(size_t n, size_t m, double eps);

/** Verify a seeded extractor by empirical testing.
 *  Draws random sources with known min-entropy, applies extractor,
 *  and tests output uniformity via chi-squared test.
 */
bool ext_seeded_verify(const SeededExt *e, size_t trials);

/*──────────────────────────────────────────────────────────────────────
 * L6: Two-Source Extractors
 *──────────────────────────────────────────────────────────────────────*/

/** Chor-Goldreich inner-product two-source extractor (1988).
 *
 *  2Ext(x,y) = Σ_i x_i·y_i (mod 2)  (inner product over GF(2))
 *
 *  For n-bit independent sources each with min-entropy > n/2,
 *  outputs 1 nearly-uniform bit.
 *
 *  This implementation leverages the Vazirani XOR lemma to output
 *  multiple bits by partitioning into blocks.
 *
 *  Theorem: If X,Y independent with H∞(X), H∞(Y) > n/2,
 *           then Δ(⟨X,Y⟩, U_1) ≤ 2^{-(H∞(X)+H∞(Y)-n-1)/2}.
 */
TwoSourceExt ext_two_source_chor_goldreich(size_t n);

/** Evaluate the Chor-Goldreich extractor on two concrete sources.
 *  Returns extracted bits (m bits stored in out, LSB-aligned).
 */
bool ext_cg_eval(const TwoSourceExt *e, const bool *x, const bool *y,
                 bool *out);

/** Bourgain two-source extractor (2005).
 *
 *  For sources with min-entropy rate > ½, extracts Ω(n) bits.
 *  Based on sum-product theorem in finite fields.
 *
 *  This implementation: encode sources as elements of GF(2^n),
 *  compute 2Ext(x,y) = x·y (field multiplication), output bits.
 *
 *  Theorem (Bourgain): For constant δ>0, if H∞(X),H∞(Y) ≥ (½+δ)n,
 *    then 2Ext extracts Ω(n) bits with ε = 2^{-Ω(n)}.
 */
TwoSourceExt ext_two_source_bourgain(size_t n, double eps);

/** Evaluate Bourgain two-source extractor on concrete inputs.
 */
bool ext_bourgain_eval(const TwoSourceExt *e, const bool *x, const bool *y,
                       bool *out);

/** Maximum output length for a two-source extractor with given
 *  source min-entropy. For Chor-Goldreich: m ≤ min(H∞(X),H∞(Y)) - n/2.
 */
size_t ext_two_source_output(const TwoSourceExt *e);

/*──────────────────────────────────────────────────────────────────────
 * L8: Extractor Graphs (expander-based constructions)
 *──────────────────────────────────────────────────────────────────────*/

/** Extractor graph: a bipartite graph G = (L∪R, E) with |L|=2^n,
 *  |R|=2^m, left-degree D=2^d, such that for any set S⊆L of size
 *  ≥ 2^k, the neighbor set is almost all of R.
 *
 *  This is equivalent to a (k,ε)-extractor.
 */
typedef struct {
    size_t n;          /**< Left side: 2^n vertices */
    size_t k;          /**< Min-entropy threshold */
    size_t d;          /**< log₂(degree) */
    size_t m;          /**< Right side: 2^m vertices */
    size_t **walks;    /**< walks[i][j] = j-th neighbor label (size_t) */
    size_t degree;     /**< actual degree = 2^d */
    double eps;        /**< Expansion error */
} ExtractorGraph;

/** Construct an extractor graph from parameters.
 *  Uses a Ramanujan expander-based construction.
 */
ExtractorGraph eg_construct(size_t n, size_t k, size_t d);

/** Verify the extractor property of the graph.
 *  For random subsets of size ≥ 2^k, check neighbor expansion.
 *  Uses spectral expansion bound for verification.
 */
bool eg_verify_expansion(const ExtractorGraph *eg);

/** Get the neighbor set size for a given left-set size.
 *  Returns fraction of right side covered.
 */
double eg_neighbor_fraction(const ExtractorGraph *eg, size_t set_size);

/** Free extractor graph resources.
 */
void eg_free(ExtractorGraph *eg);

/*──────────────────────────────────────────────────────────────────────
 * L5: Pseudorandom Generator from Extractor
 *──────────────────────────────────────────────────────────────────────*/

/** Nisan-Wigderson style PRG from a hard predicate.
 *  Given a function f: {0,1}^t → {0,1} that is hard for circuits of
 *  size S, and an (r, D, m, t)-design, outputs m pseudorandom bits.
 *
 *  NW(x)_i = f(x|_{S_i})  where S_i are the design sets.
 */
bool ext_nw_generator(const bool *seed, size_t seed_len,
                      bool (*predicate)(const bool*, size_t),
                      const size_t **design, size_t m, size_t t,
                      bool *output);

/*──────────────────────────────────────────────────────────────────────
 * L6: Block-source extractor
 *──────────────────────────────────────────────────────────────────────*/

/** Block-source extractor: handles sources structured as blocks
 *  where each block has min-entropy conditioned on previous blocks.
 *
 *  Theorem (NZ): If X = (X₁,...,X_b) is a block source with
 *  H∞(X_i | X₁,...,X_{i-1}) ≥ k per block, then can extract ≈ b·k bits.
 */
SeededExt ext_block_source(size_t n, size_t b, size_t k_per_block);

#endif /* EXT_CONSTR_H */
