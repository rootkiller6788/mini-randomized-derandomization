/*
 * list_decode_apps.h — Applications of List-Decodable Codes
 *
 * List-decodable codes have deep applications in computational
 * complexity theory and cryptography:
 *
 * L7 Applications:
 *   - Hardness amplification (Yao's XOR lemma via codes)
 *   - Randomness extractors from list-decodable codes
 *   - Soft-decision decoding (Koetter-Vardy framework)
 *   - Derandomization of BPP (via hard functions)
 *
 * L8 Advanced Topics:
 *   - Capacity-achieving list decoding
 *   - List recoverability and its uses in randomness
 *   - Concatenated codes for explicit constructions
 *
 * References:
 *   - Sudan, Trevisan, Vadhan, "Pseudorandom Generators without XOR" (JCSS 2001)
 *   - Guruswami, Umans, Vadhan, "Unbalanced Expanders and Randomness Extractors
 *     from Parvaresh-Vardy Codes" (JACM 2009)
 *   - Trevisan, "List-Decoding Using the XOR Lemma" (FOCS 2003)
 */
#ifndef LD_APPS_H
#define LD_APPS_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

/* Forward declaration for CodeParams */
struct _ld_code_params;
typedef struct _ld_code_params CodeParamsLD;
struct _ld_code_params {
    size_t n, k, d, q;
};

/* Forward declaration for Codeword */
struct _ld_codeword;
typedef struct _ld_codeword CodewordLD;
struct _ld_codeword {
    size_t *symbols;
    size_t  len;
    size_t  alphabet;
};

/* ------------------------------------------------------------------ */
/*  L4 — Fundamental Bounds (from list_decode_bounds.c)                */
/* ------------------------------------------------------------------ */

/** Hamming (sphere-packing) bound: log_q(M) ≤ n - log_q(vol(B_t)). */
double ld_hamming_bound(size_t n, size_t d, size_t q);

/** Singleton bound: k ≤ n - d + 1. */
size_t ld_singleton_bound_k(size_t n, size_t d, size_t q);

/** Singleton bound as rate: R ≤ 1 - δ + 1/n. */
double ld_singleton_bound_rate(size_t n, size_t d);

/** Plotkin bound (low-rate regime): for δ > (q-1)/q. */
double ld_plotkin_bound(size_t n, size_t d, size_t q);

/** Gilbert-Varshamov existence bound: R ≥ 1 - H_q(δ). */
double ld_gilbert_varshamov_bound(size_t n, size_t d, size_t q);

/** MRRW (LP) bound for binary codes. */
double ld_mrrw_bound(double delta);

/** Johnson bound with list-size computation (extended version). */
size_t ld_johnson_bound_extended(size_t n, size_t d, size_t q, size_t e);

/** Relaxed Singleton bound for list decoding with list size L. */
double ld_relaxed_singleton_bound(size_t n, size_t d, size_t q, size_t L);

/** Elias-Bassalygo upper bound. */
double ld_elias_bassalygo_bound(size_t n, size_t d, size_t q);

/** Capacity bound for list decoding: δ_cap = 1 - R - ε. */
double ld_list_capacity_bound(size_t n, size_t k, size_t q, double epsilon);

/** List-recovery capacity bound. */
double ld_list_recovery_capacity(size_t n, size_t k, size_t q,
                                  size_t ell, double epsilon);

/** Mutual information for q-ary symmetric channel. */
double ld_soft_capacity_mutual_info(double crossover_prob, size_t q);

/** Exact A(n,d)_q for small n. */
size_t ld_max_code_size_small_n(size_t n, size_t d, size_t q);

/** Verify all bounds for given code parameters. */
bool ld_verify_all_bounds(size_t n, size_t M, size_t d, size_t q);

/** Enumerate all codewords of a linear code (from list_decode_core.c). */
size_t ld_enumerate_codewords(const size_t *G, size_t n, size_t k, size_t q,
                               size_t **cw_out);

/** Capacity-achieving verification (from list_decode_apps.c). */
bool ld_verify_capacity_achieving(double rate, double epsilon,
                                   size_t q, size_t m);

/* ------------------------------------------------------------------ */
/*  L7 — Applications                                                  */
/* ------------------------------------------------------------------ */

/** Parameters for hardness amplification via list-decodable codes. */
typedef struct {
    size_t n;       /* input length of original function */
    double eps;     /* hardness parameter (agreement bound) */
    size_t L;       /* list-size bound for the code */
} HardnessAmpParams;

/** Parameters for randomness extraction from codes. */
typedef struct {
    size_t n;         /* input length */
    size_t k;         /* entropy of source */
    double delta;     /* bias bound */
} ExtractParams;

/**
 * Hardness amplification via list-decodable codes.
 *
 * Given a Boolean function f that is (1-ε)-hard (no circuit of size s
 * agrees with f on > 1-ε fraction of inputs), outputs an amplified
 * function g that is (1/2+δ)-hard for an appropriately chosen δ.
 *
 * Theorem (Impagliazzo 1995, Sudan-Trevisan-Vadhan 2001):
 *   A (1-ε)-hard function can be amplified to (1/2+δ)-hardness
 *   for any δ > 0 using list-decodable codes.
 *
 * @param fn          Original Boolean function values (length 2^n, truth-table).
 * @param n           Input length in bits.
 * @param hap         Hardness amplification parameters.
 * @param amplified   Output: amplified hardness threshold.
 * @return true on success.
 */
bool ld_hardness_amplification(const bool *fn, size_t n,
                               const HardnessAmpParams *hap,
                               double *amplified);

/**
 * Construct a randomness extractor from a list-decodable code.
 * Uses the Trevisan extractor paradigm.
 *
 * @param source   Weak random source bits (length slen).
 * @param output   Extracted near-uniform bits (length olen).
 * @return true on success.
 */
bool ld_extractor_from_codes(const CodeParamsLD *cp,
                             const bool *source, size_t slen,
                             bool *output, size_t olen);

/**
 * Soft-decision decoding: given reliability values for each symbol,
 * find the most likely codeword (MAP decoding).
 *
 * Based on the Koetter-Vardy (2003) algebraic soft-decision decoder.
 */
bool ld_soft_decoding(const CodewordLD *received,
                      const double *reliability,
                      size_t *decoded);

/* ------------------------------------------------------------------ */
/*  L4 — Fundamental Bounds (Additional)                               */
/* ------------------------------------------------------------------ */

/**
 * Capacity achievability check: for a given q, error fraction δ,
 * and target list size L, determines if list decoding is
 * combinatorially possible.
 *
 * Theorem (Zyablov-Pinsker 1981): For any ε > 0, there exist codes
 * of rate R = 1 - δ - ε that are list-decodable up to δ with
 * list size O(1/ε).
 *
 * @return Achievable rate for parameters (q, δ, L).
 */
double ld_capacity_achievability(size_t q, double delta, size_t L);

/**
 * Plotkin bound: for a code of length n and distance d over alphabet
 * of size q, if d > n·(q-1)/q then the code size is bounded by d/(d-n+q/α).
 *
 * Simplified version: |C| ≤ d / (d - n·(q-1)/q).
 */
bool ld_plotkin_bound_test(size_t n, size_t d, size_t q, size_t code_size);

/**
 * Singleton bound: for an [n,k,d]_q code, k ≤ n - d + 1.
 */
size_t ld_singleton_bound(size_t n, size_t d, size_t q);

/**
 * Test if the MDS conjecture holds for parameters (n,k,q).
 * MDS conjecture: there exists an [n,k,n-k+1]_q code (linear MDS)
 * iff n ≤ q+1 (for prime q) or n ≤ q+2 (for q=2^m with n=3 or n=q-1).
 *
 * @return true iff an MDS code is known to exist for given parameters.
 */
bool ld_mds_conjecture_verify(size_t n, size_t k, size_t q);

/* ------------------------------------------------------------------ */
/*  L8 — Advanced Applications                                         */
/* ------------------------------------------------------------------ */

/**
 * Concatenated code encoding: outer code (list-decodable) + inner code
 * (small alphabet).  Classic Forney construction.
 *
 * @param outer_n    Outer code block length.
 * @param inner_n    Inner code block length.
 * @param msg        Message to encode.
 * @param codeword   Output codeword of length outer_n * inner_n.
 * @return true on success.
 */
bool ld_concatenated_encode(const CodeParamsLD *outer,
                            const CodeParamsLD *inner,
                            const size_t *msg,
                            size_t *codeword);

/**
 * List-recoverability: given a set S_i ⊆ [q] of size ℓ for each i ∈ [n],
 * find all codewords c such that c_i ∈ S_i for at least α·n positions.
 * Essential for extractor and cryptographic constructions.
 */
size_t ld_list_recover(const CodewordLD **input_sets,
                       size_t *set_sizes,
                       size_t n, size_t ell, double alpha,
                       CodewordLD **decoded, size_t max_list);

#endif /* LD_APPS_H */
