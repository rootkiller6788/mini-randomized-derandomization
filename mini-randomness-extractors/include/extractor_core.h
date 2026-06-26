#ifndef EXT_CORE_H
#define EXT_CORE_H

/**
 * mini-randomness-extractors — Core Definitions
 *
 * L1 Definitions: Extractor, WeakSource, SeededExt, TwoSourceExt, Disperser
 * L2 Core Concepts: min-entropy, statistical distance, seed length, output length
 * L3 Mathematical Structures: probability simplex, total variation distance, entropy measures
 *
 * References:
 *   Nisan & Zuckerman (1996) "Randomness is Linear in Space"
 *   Trevisan (2001) "Extractors and Pseudorandom Generators"
 *   Shaltiel (2002) "Recent Developments in Extractors"
 *   Vadhan (2012) "Pseudorandomness" §6
 */

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>

/*──────────────────────────────────────────────────────────────────────
 * L1: Core type definitions
 *──────────────────────────────────────────────────────────────────────*/

/** A (k,ε)-extractor Ext : {0,1}^n × {0,1}^d → {0,1}^m.
 *  For every distribution X over {0,1}^n with H∞(X) ≥ k,
 *  Δ(Ext(X, U_d), U_m) ≤ ε.
 *
 *  n = source length, d = seed length, m = output length, eps = error ε
 */
typedef struct {
    size_t n;      /**< Source length (bits) */
    size_t d;      /**< Seed length (bits) */
    size_t m;      /**< Output length (bits) */
    double eps;    /**< Statistical error ε ∈ [0,1] */
} Extractor;

/** A k-source: a distribution X over {0,1}^n with H∞(X) ≥ k.
 *  min_entropy = k = -log₂(max_x Pr[X=x]).
 */
typedef struct {
    size_t n;           /**< Domain size (bits) */
    double min_entropy; /**< Guaranteed min-entropy (bits) */
} WeakSource;

/** A seeded extractor with an embedded evaluation function.
 *  n = source length, m = output length.
 *  The function pointer ext(src, seed, out, out_len) evaluates the extractor.
 *  Parameter interpretation:
 *    src     — n-bit weak source string
 *    seed    — seed index (encodes the d-bit seed as an integer for small d)
 *    out     — output buffer (m bits written)
 *    out_len — must equal m
 */
typedef struct {
    size_t n;        /**< Source length (bits) */
    size_t m;        /**< Output length (bits) */
    bool (*ext)(const bool *src, size_t seed, bool *out, size_t out_len);
} SeededExt;

/** A two-source extractor 2Ext : {0,1}^n × {0,1}^n → {0,1}^m.
 *  For independent sources X,Y each with min-entropy ≥ k,
 *  Δ(2Ext(X,Y), U_m) ≤ ε.
 */
typedef struct {
    size_t n;        /**< Each source length (bits) */
    size_t m;        /**< Output length (bits) */
    double eps;      /**< Statistical error ε */
} TwoSourceExt;

/*──────────────────────────────────────────────────────────────────────
 * L4: Fundamental bounds (see extractor_core.c for derivations)
 *──────────────────────────────────────────────────────────────────────*/

/** Create an extractor with given parameters.
 *  Validates: d ≤ n, m ≤ k + 2·log₂(eps) + O(1), eps ∈ (0,1].
 *  Theorem (RTD, NZ): m ≤ k + 2·log₂(ε) + 2, d ≥ log(n-k) + 2·log(1/ε) + O(1).
 */
Extractor ext_create(size_t n, size_t d, size_t m, double eps);

/** Evaluate the extractor on a concrete weak source and seed.
 *  source_bits: n-bit string drawn from a k-source
 *  seed:        d-bit uniformly random seed
 *  out:         output buffer of m bits
 *  Returns true iff extraction succeeded within error bound.
 */
bool ext_eval(const Extractor *e, const bool *source_bits,
              const bool *seed, bool *out);

/** Compute the theoretical error bound ε for extractor with parameters
 *  (n, d, m) applied to a source of min-entropy k.
 *  LHL bound: ε ≤ 2^{-(k-m)/2}.
 */
double ext_error_bound(const Extractor *e, size_t n, double k);

/** Minimum seed length required to achieve error ε for source length n
 *  with min-entropy k.
 *  d ≥ max(1, log₂(n-k) + 2·log₂(1/ε) + O(1)).
 */
size_t ext_seed_len(size_t n, double min_entropy, double eps);

/** Maximum output length achievable from source length n with
 *  min-entropy k (before hitting ε = 2^{-k/2} threshold).
 *  m ≤ k - 2·log₂(1/ε) - O(1).
 */
size_t ext_output_len(size_t n, double min_entropy);

/** Validate extractor parameters for consistency.
 *  Checks: n>0, d≤n, m>0, eps∈(0,1], m≤n.
 */
bool ext_is_valid(const Extractor *e);

/** Free extractor resources (no-op for parameter-only extractors). */
void ext_free(Extractor *e);

/*──────────────────────────────────────────────────────────────────────
 * L2: Entropy measures
 *──────────────────────────────────────────────────────────────────────*/

/** Compute min-entropy H∞(X) from a sequence of samples.
 *  H∞(X) = -log₂(max_x Pr[X=x]).
 *  dist: array of nsamples * n bits (each sample is n bits packed)
 *  nsamples: number of samples
 *  n: bit-length of each sample
 *  Returns H∞ in bits, or -1.0 on error.
 */
double ext_min_entropy(const bool *dist, size_t nsamples, size_t n);

/** Compute min-entropy from explicit frequency table.
 *  freq[i] = count of samples equal to value i.
 *  nvals: number of possible values (= 2^n for n-bit samples)
 *  nsamples: total number of samples (sum of freq)
 *  Returns H∞ in bits.
 */
double ext_min_entropy_from_freq(const size_t *freq, size_t nvals,
                                  size_t nsamples);

/** Compute Shannon entropy H(X) = -Σ Pr[X=x] · log₂ Pr[X=x].
 */
double ext_shannon_entropy(const size_t *freq, size_t nvals,
                            size_t nsamples);

/** Compute Rényi entropy of order α.
 *  H_α(X) = (1/(1-α)) · log₂(Σ_x Pr[X=x]^α).
 *  α=2 gives collision entropy; α→∞ gives min-entropy.
 */
double ext_renyi_entropy(const size_t *freq, size_t nvals,
                          size_t nsamples, double alpha);

/** Compute statistical (total variation) distance between two
 *  distributions over the same domain.
 *  Δ(P,Q) = ½ Σ_x |P(x) - Q(x)| = max_S |P(S) - Q(S)|.
 *  freq1, freq2: frequency tables over the same nvals values
 *  total1, total2: sample counts
 *  Returns Δ ∈ [0,1].
 */
double ext_statistical_distance(const size_t *freq1, const size_t *freq2,
                                 size_t nvals, size_t total1, size_t total2);

/** Compute Hellinger distance: H(P,Q) = (1/√2) · √(Σ_x (√P(x) - √Q(x))²).
 *  Bounded by: H² ≤ Δ ≤ √2·H.
 */
double ext_hellinger_distance(const size_t *freq1, const size_t *freq2,
                               size_t nvals, size_t total1, size_t total2);

/** Compute Kullback-Leibler divergence D_KL(P||Q).
 *  D_KL(P||Q) = Σ_x P(x)·log(P(x)/Q(x)).
 *  Not symmetric; returns ∞ if Q(x)=0 and P(x)>0.
 */
double ext_kl_divergence(const size_t *freq1, const size_t *freq2,
                          size_t nvals, size_t total1, size_t total2);

/*──────────────────────────────────────────────────────────────────────
 * Utility: bit-vector helpers
 *──────────────────────────────────────────────────────────────────────*/

/** Pack n booleans into a bit-vector (LSB first).
 *  bits: array of n bools, vec: output size_t array (ceil(n/wordsize) words).
 */
void ext_bits_pack(const bool *bits, size_t n, size_t *vec);

/** Unpack a bit-vector into boolean array.
 */
void ext_bits_unpack(const size_t *vec, size_t n, bool *bits);

/** Count number of set bits in a bool array of length n.
 */
size_t ext_popcount_bool(const bool *bits, size_t n);

/** Hamming weight of a size_t value.
 */
size_t ext_popcount(size_t x);

/** LCG PRNG step (exposed for cross-module use).
 *  Returns next pseudorandom 64-bit value and updates state.
 */
uint64_t ext_lcg_next(uint64_t *state);

/** LCG PRNG: uniform double in [0,1). */
double ext_lcg_double(uint64_t *state);

/** Generate a uniformly random bool array of length n.
 *  Uses a simple LCG for reproducibility; for cryptographic use,
 *  replace with system entropy source.
 */
void ext_random_bits(bool *bits, size_t n, uint64_t *seed);

/** Fill a frequency table by sampling from a biased source.
 *  bias: probability of 1 for each independent bit.
 *  nsamples: number of n-bit samples to draw.
 *  freq: output frequency table of size 2^n (WARNING: n must be ≤ 16).
 */
void ext_sample_distribution(double bias, size_t n, size_t nsamples,
                              size_t *freq, uint64_t *seed);

/** Fill a frequency table for a "block source" where each sample
 *  consists of n independent bits, but with bit-wise varying biases.
 *  biases: array of n bias values ∈ [0,1].
 */
void ext_sample_block_source(const double *biases, size_t n,
                              size_t nsamples, size_t *freq, uint64_t *seed);

/** Compute the entropy deficit (gap to uniform) in bits.
 *  deficit = n - H∞(X).  Uniform distribution → deficit = 0.
 */
double ext_entropy_deficit(size_t n, double min_entropy);

/** Verify that an empirical distribution's min-entropy is ≥ k
 *  with confidence 1-δ using Chernoff bounds.
 */
bool ext_verify_min_entropy(const size_t *freq, size_t nvals,
                             size_t nsamples, double k, double delta);

/* GF(2) field operations (L3) */
void ext_gf2_add(const bool *a, const bool *b, size_t n, bool *result);
bool ext_gf2_mul(const bool *a, const bool *b, size_t n, bool *result);

/* Bit conversion utilities */
size_t ext_bits_to_size(const bool *bits, size_t n);
void ext_size_to_bits(size_t val, size_t n, bool *bits);
void ext_bits_print(const bool *bits, size_t n, const char *label);
void ext_bits_copy(const bool *src, bool *dst, size_t n);
void ext_bits_xor(const bool *a, const bool *b, bool *result, size_t n);
bool ext_bits_eq(const bool *a, const bool *b, size_t n);

#endif /* EXT_CORE_H */
