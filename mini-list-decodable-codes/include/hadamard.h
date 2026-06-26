/*
 * hadamard.h — Hadamard Codes and the Goldreich-Levin Algorithm
 *
 * The Hadamard code encodes an n-bit message x ∈ {0,1}^n as the
 * 2^n-bit truth table of the linear function y ↦ ⟨x, y⟩ mod 2.
 *
 * Parameters: n = 2^k, k = message length, d = n/2, rate = k/2^k.
 * This is a [2^k, k, 2^{k-1}]_2 code — extremely low rate, but:
 *   - Locally decodable with 2 queries
 *   - List-decodable to radius 1/2 - ε (Goldreich-Levin theorem)
 *
 * The Goldreich-Levin (GL) algorithm is a cornerstone result:
 * given oracle access to f: {0,1}^k → {0,1}, it finds all x
 * such that f(y) = ⟨x,y⟩ for at least 1/2 + ε fraction of y.
 * This is list-decoding the Hadamard code.
 *
 * Applications: hard-core predicate extraction, learning parity
 * with noise (LPN), cryptography, property testing.
 *
 * References:
 *   - Goldreich & Levin, STOC 1989
 *   - Goldreich, Rubinfeld, Sudan, "Learning Polynomials with Queries:
 *     The Highly Noisy Case" (FOCS 1995)
 *   - Kushilevitz & Mansour, "Learning Decision Trees Using the
 *     Fourier Spectrum" (STOC 1991)
 */
#ifndef HADAMARD_H
#define HADAMARD_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

/* ------------------------------------------------------------------ */
/*  L1 — Hadamard Code Structure                                       */
/* ------------------------------------------------------------------ */

/** Hadamard code H_k: maps {0,1}^k → {0,1}^{2^k}. */
typedef struct {
    size_t k;        /* message length */
    size_t n;        /* block length = 2^k */
    size_t d;        /* minimum distance = 2^{k-1} */
} HadamardCode;

/* ------------------------------------------------------------------ */
/*  L2 — Hadamard Code Operations                                      */
/* ------------------------------------------------------------------ */

/** Create a Hadamard code of dimension k. */
HadamardCode *had_create(size_t k);

/** Encode: given k-bit message x, output n-bit codeword where
 *  codeword[y] = ⟨x, y⟩ mod 2 for all y ∈ {0,1}^k. */
void had_encode(const HadamardCode *hc, const bool *msg, bool *codeword);

/** Free a Hadamard code. */
void had_free(HadamardCode *hc);

/* ------------------------------------------------------------------ */
/*  L5 — Goldreich-Levin List-Decoding Algorithm                       */
/* ------------------------------------------------------------------ */

/**
 * Goldreich-Levin algorithm: find all x ∈ {0,1}^k such that
 * oracle(y) = ⟨x, y⟩ for at least (1/2 + ε) fraction of y ∈ {0,1}^k.
 *
 * Theorem (Goldreich-Levin 1989):
 *   There exists a probabilistic algorithm that, given oracle access
 *   to f: {0,1}^k → {0,1}, runs in time poly(k, 1/ε) and outputs
 *   a list of size O(1/ε²) containing all x with agreement ≥ 1/2+ε.
 *
 * @param oracle      Function f: {0,1}^k → {0,1}.
 * @param epsilon     Agreement advantage over random guessing.
 * @param decoded_out Output: list of candidate messages (ListSize × k).
 * @param max_list    Maximum candidates to return.
 * @return Number of candidate messages found.
 */
size_t had_goldreich_levin(const HadamardCode *hc,
                            bool (*oracle)(const bool *y),
                            double epsilon,
                            bool **decoded_out,
                            size_t max_list);

/**
 * Deterministic variant of Goldreich-Levin using explicit
 * ε-biased sets (Naor-Naor 1993) instead of random sampling.
 * Runs in time poly(k, 1/ε) deterministically.
 */
size_t had_goldreich_levin_deterministic(const HadamardCode *hc,
                                          bool (*oracle)(const bool *y),
                                          double epsilon,
                                          bool **decoded_out,
                                          size_t max_list);

/* ------------------------------------------------------------------ */
/*  L4 — Bounds for Hadamard Codes                                     */
/* ------------------------------------------------------------------ */

/**
 * List-decoding radius of Hadamard code: δ = 1/2 - ε.
 * Any word has at most 1/(2ε)² codewords within radius 1/2 - ε.
 */
double had_list_decoding_radius(double epsilon);

/**
 * List-size bound (Johnson bound for Hadamard):
 *   L_J ≤ 1 / (1 - 2δ)² = 1 / (4ε²).
 */
size_t had_list_size_bound(double epsilon);

/* ------------------------------------------------------------------ */
/*  L7 — Application: Hard-Core Predicate Extraction                   */
/* ------------------------------------------------------------------ */

/**
 * Given a one-way function f, extract a hard-core predicate using
 * the Goldreich-Levin theorem.
 *
 * If f is a one-way function, then g(x,r) = (f(x), r) has h(x,r) = ⟨x,r⟩
 * as a hard-core predicate.
 *
 * @param f_input   Input to the one-way function.
 * @param r         Random bits (length k).
 * @return ⟨f_input, r⟩ mod 2 — the hard-core bit.
 */
bool had_hard_core_bit(const bool *f_input, const bool *r, size_t k);

/* ------------------------------------------------------------------ */
/*  L8 — Fourier Analysis over Hadamard                                */
/* ------------------------------------------------------------------ */

/**
 * Compute the Fourier coefficient of a Boolean function at point α:
 *   f̂(α) = (1/2^k) · Σ_{y∈{0,1}^k} (-1)^{f(y) ⊕ ⟨α,y⟩}.
 *
 * Large Fourier coefficients correspond to good linear approximations
 * — precisely what Goldreich-Levin finds.
 */
double had_fourier_coefficient(bool (*f)(const bool *y), size_t k,
                                const bool *alpha);

/**
 * Find all α with |f̂(α)| ≥ τ (the "heavy Fourier coefficients").
 * Equivalent to list-decoding Hadamard with agreement 1/2 + τ/2.
 */
size_t had_heavy_fourier_coefficients(bool (*f)(const bool *y),
                                       size_t k, double tau,
                                       bool **alphas_out, size_t max_list);

#endif /* HADAMARD_H */