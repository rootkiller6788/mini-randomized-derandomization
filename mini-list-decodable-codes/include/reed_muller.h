/*
 * reed_muller.h — Reed-Muller Codes
 *
 * Reed-Muller (RM) codes are a classical family of binary codes
 * based on evaluations of multivariate polynomials over GF(2).
 *
 * For parameters r (degree) and m (number of variables), the code
 * RM(r, m) has:
 *   n = 2^m         (block length)
 *   k = Σ_{i=0}^{r} C(m, i)  (dimension)
 *   d = 2^{m-r}     (minimum distance)
 *
 * RM codes are locally testable, locally decodable, and form the
 * basis for many complexity-theoretic constructions including
 * the PCP theorem and hardness amplification.
 *
 * References:
 *   - Reed, IRE Trans. Info. Theory 4(4):38-49, 1954
 *   - Muller, IRE Trans. EC-3(3):6-12, 1954
 *   - Gopalan, Klivans, Zuckerman, "List-Decoding Reed-Muller Codes
 *     over Small Fields" (STOC 2008)
 *   - Kaufman, Lovett, Porat, "Weight Distribution and List-Decoding
 *     Size of Reed-Muller Codes" (IEEE Trans. IT 2012)
 */
#ifndef REED_MULLER_H
#define REED_MULLER_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

/* ------------------------------------------------------------------ */
/*  L1 — Reed-Muller Code Structure                                    */
/* ------------------------------------------------------------------ */

/** Reed-Muller code RM(r,m) over GF(2). */
typedef struct {
    size_t r;        /* maximum degree of evaluated polynomials */
    size_t m;        /* number of variables */
    size_t n;        /* block length = 2^m */
    size_t k;        /* dimension = Σ_{i=0}^r C(m,i) */
    size_t d;        /* minimum distance = 2^{m-r} */
    size_t *truth_table; /* truth table of encoder (optional, size n·k) */
} ReedMullerCode;

/* ------------------------------------------------------------------ */
/*  L2 — RM Code Operations                                            */
/* ------------------------------------------------------------------ */

/** Create a Reed-Muller code RM(r,m).
 *  @return NULL if r > m (invalid parameters). */
ReedMullerCode *rm_create(size_t r, size_t m);

/** Encode a message into an RM codeword.
 *  msg has k bits (coefficients of the multilinear polynomial).
 *  codeword_out has n = 2^m bits (evaluations at all points in {0,1}^m). */
void rm_encode(const ReedMullerCode *rm, const bool *msg,
               bool *codeword_out);

/** Free an RM code. */
void rm_free(ReedMullerCode *rm);

/** Compute binomial coefficient C(n,k). */
size_t rm_binomial(size_t n, size_t k);

/** Compute dimension k = Σ C(m,i) for i=0..r. */
size_t rm_dimension(size_t r, size_t m);

/* ------------------------------------------------------------------ */
/*  L5 — List-Decoding Algorithms for RM Codes                         */
/* ------------------------------------------------------------------ */

/**
 * List-decode an RM(r,m) code using the Gopalan-Klivans-Zuckerman
 * (GKZ 2008) algorithm.
 *
 * The algorithm works by reducing multivariate list-decoding to
 * univariate list-decoding along random lines.
 *
 * @param received      Received word (length n = 2^m, bits).
 * @param error_frac    Error fraction δ.
 * @param decoded_out   Output: array of candidate messages (L × k bits).
 * @param max_list      Maximum list size.
 * @return Number of candidates found.
 */
size_t rm_list_decode_gkz(const ReedMullerCode *rm,
                           const bool *received,
                           double error_frac,
                           bool **decoded_out,
                           size_t max_list);

/**
 * Local list-decoding of RM codes.
 * Given oracle access to a received word, find all polynomials
 * that agree on at least (1-δ) fraction of points, querying
 * only a small number of coordinates.
 *
 * @param oracle        Function that gives value at point x ∈ {0,1}^m.
 * @param decoded_out   Output candidate messages.
 * @param max_list      Max list size.
 * @return Number of candidates.
 */
size_t rm_local_list_decode(const ReedMullerCode *rm,
                             bool (*oracle)(const bool *x),
                             double error_frac,
                             bool **decoded_out,
                             size_t max_list);

/* ------------------------------------------------------------------ */
/*  L4 — List-Decoding Bounds for RM Codes                             */
/* ------------------------------------------------------------------ */

/**
 * List-decoding radius for RM(r,m) codes (GKZ bound).
 * For δ < 1 - sqrt{r/(2^{m-1})}, list size is polynomially bounded.
 */
double rm_list_decoding_radius(const ReedMullerCode *rm);

/**
 * List-size bound for RM(r,m) codes.
 * Theorem (GKZ): For δ_opt = 2·(2^{-r}), list size ≤ O(n²).
 */
size_t rm_list_size_bound(const ReedMullerCode *rm, double error_frac);

/**
 * Verify that a word can be uniquely decoded (within d/2 errors). */
bool rm_is_uniquely_decodable(const ReedMullerCode *rm, size_t errors);

/* ------------------------------------------------------------------ */
/*  L6 — Canonical Problem: RM Decoding                                */
/* ------------------------------------------------------------------ */

/**
 * Majority-logic decoder for RM codes (Reed's algorithm).
 * This is a simple unique-decoding algorithm up to d/2 errors.
 * Acts as a baseline for comparison with list-decoding.
 *
 * @return true on success.
 */
bool rm_majority_logic_decode(const ReedMullerCode *rm,
                               const bool *received,
                               bool *decoded_out);

#endif /* REED_MULLER_H */