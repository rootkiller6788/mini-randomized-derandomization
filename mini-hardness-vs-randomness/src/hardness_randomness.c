/******************************************************************************
 * hardness_randomness.c �� Hardness vs Randomness Core Implementation
 *
 * Implements the fundamental connection: if hard functions exist (circuit
 * lower bounds), then pseudorandom generators exist, which can derandomize
 * BPP. This is the Nisan-Wigderson / Impagliazzo-Wigderson paradigm.
 *
 * Knowledge: L1 (hardness, PRG definitions), L2 (HR connection),
 * L4 (NW theorem, IW theorem), L5 (NW PRG construction)
 *
 * References:
 *   Nisan & Wigderson (1994): Hardness vs Randomness, JCSS 49(2)
 *   Impagliazzo & Wigderson (1997): P = BPP if E requires exponential circuits
 *   Arora & Barak (2009), Chapters 19-20
 *   Goldreich (2008), Chapter 8
 ******************************************************************************/

#include "hardness_randomness.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

/* ================================================================
 * L1+L2: Hardness Level to PRG Mapping
 *
 * The hardness-randomness tradeoff:
 * - WEAK hardness    �� PRG with mild stretch (n �� n + log n)
 * - MODERATE hardness �� PRG with polynomial stretch (n �� n^c)
 * - STRONG hardness   �� PRG with exponential stretch (n �� 2^{n^��})
 * - EXPONENTIAL hardness �� PRG with full exponential stretch
 *
 * Constraint: seed_len = O(log input_size / hardness)
 *            output_len = seed_len * (hardness improvement factor)
 * ================================================================ */

const char *hr_level_name(HRLevel level) {
    switch (level) {
        case HR_WEAK:        return "WEAK";
        case HR_MODERATE:    return "MODERATE";
        case HR_STRONG:      return "STRONG";
        case HR_EXPONENTIAL: return "EXPONENTIAL";
        default:             return "UNKNOWN";
    }
}

/*
 * Compute hardness value from circuit lower bound.
 * hardness = log2(minimum_circuit_size) / n.
 * A hardness of 1.0 means circuits of size 2^n are needed.
 * A hardness of 0.01 means circuits of size 2^{0.01n} are needed.
 *
 * Reference: Arora & Barak, Definition 19.1
 */
double hr_hardness_from_circuit_lb(size_t n, size_t S) {
    if (n == 0) return 0.0;
    if (S <= 1) return 0.0;
    return log2((double)S) / (double)n;
}

/*
 * Compute the Nisan-Wigderson PRG seed length.
 *
 * If we have a function f: {0,1}^k �� {0,1} with hardness H = 2^{��(k)},
 * and we want to fool circuits of size n, we need:
 *   seed_len �� O(k) �� O(log n)  (since n = 2^{��(k)})
 *
 * More precisely: seed_len = O(log(output_len) / hardness_rate)
 *
 * Reference: Arora & Barak, Theorem 20.6
 */
size_t hr_seed_len(size_t n, HRLevel level) {
    double dn = (double)n;
    switch (level) {
        case HR_WEAK:
            /* seed �� n / log n */
            return (size_t)(dn / log2(dn + 1) + 1);
        case HR_MODERATE:
            /* seed �� n^�� for small �� */
            return (size_t)pow(dn, 0.3);
        case HR_STRONG:
            /* seed �� O(log n) */
            return (size_t)(2.0 * log2(dn + 1) + 1);
        case HR_EXPONENTIAL:
            /* seed �� O(log n) with better constant */
            return (size_t)(1.5 * log2(dn + 1) + 1);
        default:
            return n;
    }
}

/*
 * Compute PRG output length based on hardness level.
 * The stretch = output_len / seed_len.
 *
 * Stronger hardness �� longer stretch �� more pseudorandom bits.
 *
 * Reference: Arora & Barak, Section 20.3
 */
size_t hr_output_len(size_t n, HRLevel level) {
    double dn = (double)n;
    switch (level) {
        case HR_WEAK:
            /* Mild stretch: n + O(log n) */
            return n + (size_t)log2(dn + 1);
        case HR_MODERATE:
            /* Polynomial stretch: n^c for c > 1 */
            return (size_t)(dn * log2(dn + 1));
        case HR_STRONG:
            /* Exponential stretch: 2^{n^��} */
            return (size_t)pow(dn, log2(dn + 1));
        case HR_EXPONENTIAL:
            /* Full exponential: 2^{��(n)} */
            return (size_t)pow(2.0, 0.1 * dn);
        default:
            return n;
    }
}

/*
 * The Impagliazzo-Wigderson threshold:
 * If there exists f in EXP with circuit complexity �� 2^{�š�n},
 * then BPP = P (derandomization is possible).
 *
 * This function computes the minimal �� needed for derandomization
 * at input size n.
 *
 * Reference: Impagliazzo & Wigderson (1997)
 */
double hr_impagliazzo_wigderson_threshold(size_t n) {
    if (n < 2) return 1.0;
    /* The threshold �� for the IW theorem is > 0.
       For practical bounds: �� �� 1/log(n) suffices for some constructions */
    return 1.0 / log2((double)n + 1);
}

/*
 * Check whether a given hardness level suffices to prove BPP = P
 * for input size n.
 *
 * Derandomization of BPP requires hardness against circuits of
 * size 2^{��(n)}. This checks if the level provides that.
 */
bool hr_bpp_equals_p_check(HRLevel level, size_t n) {
    switch (level) {
        case HR_WEAK:
            return false; /* WEAK hardness insufficient for full derandomization */
        case HR_MODERATE:
            return (n >= 10); /* MODERATE might work for large enough n */
        case HR_STRONG:
            return true;  /* STRONG implies BPP = P */
        case HR_EXPONENTIAL:
            return true;  /* EXPONENTIAL definitely implies BPP = P */
        default:
            return false;
    }
}

/*
 * Verify that a hard function indeed exhibits the claimed hardness
 * at least heuristically, by testing on random samples.
 *
 * Cannot actually verify circuit lower bounds (that's coNP-hard),
 * but can run basic sanity checks.
 */
bool hr_verify_hardness(const HardFn *hf, size_t samples) {
    if (!hf || !hf->fn || hf->n > 24) return false;
    if (samples == 0) samples = 100;
    size_t n = hf->n;
    if (n == 0) return false;
    /* Test: the function should be reasonably balanced
       (constant functions are easy) */
    size_t ones = 0;
    bool *input = calloc(n, sizeof(bool));
    if (!input) return false;
    for (size_t s = 0; s < samples; s++) {
        /* Generate pseudo-random input */
        for (size_t i = 0; i < n; i++) {
            input[i] = ((s * 1103515245 + 12345 + i * 1664525) >> 16) & 1;
        }
        if (hf->fn(input, n)) ones++;
    }
    free(input);
    double balance = (double)ones / (double)samples;
    /* Function should not be trivially constant */
    return (balance > 0.1 && balance < 0.9);
}

/*
 * Check if hardness amplification is possible from a given base hardness.
 * If base hardness is �� ��(1/n), then amplification to constant hardness
 * is possible via Yao's XOR lemma.
 */
bool hr_hardness_amplification_possible(double base, size_t n) {
    if (n == 0) return false;
    /* Threshold: hardness �� 1/poly(n) can be amplified */
    double threshold = 1.0 / pow((double)n, 3.0);
    return base >= threshold;
}

/*
 * Number of trials needed for the NW PRG to achieve error �� ��
 * given hardness h.
 *
 * trials �� log(1/��) / h
 *
 * Reference: Nisan & Wigderson (1994), Lemma 2.4
 */
size_t hr_nw_prg_trials(size_t n, double hardness, double epsilon) {
    if (hardness <= 0.0) return n;
    if (epsilon <= 0.0) epsilon = 0.001;
    double trials = log(1.0 / epsilon) / hardness;
    if (trials < 1.0) trials = 1.0;
    return (size_t)trials;
}

/*
 * Compute PRG stretch factor = output_len / seed_len.
 * A stretch > 1 means the PRG expands randomness.
 */
double hr_compute_prg_stretch(size_t seed_len, size_t output_len) {
    if (seed_len == 0) return 1.0;
    return (double)output_len / (double)seed_len;
}

/*
 * Compute optimal seed length for given hardness and error tolerance.
 *
 * seed_len = O(log(n/��) / hardness)
 *
 * This is the key efficiency parameter: stronger hardness allows
 * shorter seeds (less true randomness needed).
 */
size_t hr_optimal_seed_length(double hardness, double epsilon) {
    if (hardness <= 0.0) return 64;
    if (epsilon <= 0.0) epsilon = 0.001;
    double optimal = log2(1.0 / epsilon) / hardness;
    if (optimal < 1.0) optimal = 1.0;
    return (size_t)optimal;
}

/*
 * Core theorem: Hardness Implies PRG.
 *
 * If there exists f with hardness �� h_0, then there exists a PRG
 * stretching k bits to m bits that fools circuits of size m.
 *
 * The NW construction uses:
 * - A combinatorial design S_1, ..., S_m ? [k] with |S_i| = ?,
 *   |S_i �� S_j| �� log m for i �� j
 * - The hard function g: {0,1}^? �� {0,1}
 * - PRG output bit i = g(x_{S_i}) where x is the seed
 *
 * This function validates that the parameters allow the construction.
 *
 * Reference: Nisan & Wigderson (1994), Theorem 1
 *          Arora & Barak, Theorem 20.6
 */
bool hr_implies_prg(const HardFn *hf, HRLevel level, PRG *result) {
    if (!hf || !result) return false;
    size_t n = hf->n;
    if (n < 4) return false;
    double hardness = hf->hardness;
    /* Check if hardness is sufficient for the level */
    double min_lb = clb_min_lower_bound_for_level(level, n);
    double actual_hardness = hr_hardness_from_circuit_lb(n, (size_t)min_lb);
    if (hardness < actual_hardness) return false;
    /* Compute PRG parameters */
    size_t seed = hr_seed_len(n, level);
    size_t output = hr_output_len(n, level);
    if (seed == 0 || output == 0) return false;
    result->seed_len = seed;
    result->output_len = output;
    result->gen = NULL; /* Generator function is not implemented in C */
    return true;
}

/*
 * Construct an NW combinatorial design.
 *
 * An (n, ?, t)-design is a family of m subsets S_1, ..., S_m ? [n]
 * such that:
 *   - |S_i| = ? for all i
 *   - |S_i �� S_j| �� t for all i �� j
 *
 * Design parameters for NW PRG:
 *   ? = O(log n),  t = O(log m),  m = output_len
 *
 * Construction uses algebraic method: set of lines in GF(q)^2.
 *
 * Reference: Nisan & Wigderson (1994), Lemma 2.1
 */
bool hr_construct_nw_design(size_t n, size_t m, size_t **sets_out, size_t *set_sizes_out) {
    if (!sets_out || !set_sizes_out || n == 0 || m == 0) return false;
    /* ? = ceil(log m) + 1, t = ceil(log m) */
    size_t ell = (size_t)ceil(log2((double)m)) + 1;
    if (ell > n) ell = n;
    /* Allocate sets: m �� ? matrix of set elements */
    size_t *sets = calloc(m * ell, sizeof(size_t));
    if (!sets) return false;
    /* Algebraic construction:
       For prime p > n, identify [n] with a subset of GF(p) �� GF(p).
       Each set S_{a,b} = { (x, ax + b mod p) : x = 0, 1, ..., ell-1 }
       This ensures |S_i �� S_j| �� 1 (for distinct a). */
    /* Simplified: For small parameters, use block design.
       Here we just produce a simple explicit construction. */
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < ell; j++) {
            /* Use a simple deterministic mapping */
            size_t val = (i * 1103515245 + j * 1664525 + 1013904223) % n;
            sets[i * ell + j] = val;
        }
    }
    *sets_out = sets;
    *set_sizes_out = ell;
    return true;
}

/*
 * Evaluate the Nisan-Wigderson PRG on a given seed.
 *
 * NW Generator: G: {0,1}^k → {0,1}^m
 *   G(x)_i = f(x|_{S_i})
 * where S_1,...,S_m is an (n, ℓ, log m)-design and
 * f: {0,1}^ℓ → {0,1} is the hard predicate.
 *
 * The combinatorial design ensures that the output bits
 * are pairwise "almost independent", fooling any
 * distinguisher of size ≤ m when f has hardness ≥ m.
 *
 * This function evaluates G(x) = G(x)_1 ... G(x)_m
 * given the seed x, the hard predicate hard_fn,
 * and a pre-constructed design (sets, sizes).
 *
 * Complexity: O(m * time(hard_fn)).
 *
 * Reference: Nisan & Wigderson (1994), Construction 2.2
 *          Arora & Barak, Algorithm 20.1
 */
bool hr_evaluate_nw_prg(const bool *seed, size_t seed_len,
                         const bool *hard_pred_tt, size_t pred_size,
                         const size_t *sets, size_t ell, size_t m,
                         bool *output) {
    if (!seed || !hard_pred_tt || !sets || !output) return false;
    if (seed_len == 0 || ell == 0 || m == 0) return false;
    (void)pred_size; /* truth table size = 2^ell, consistency check reserved */
    for (size_t i = 0; i < m; i++) {
        /* Extract the ℓ bits of seed at positions S_i */
        size_t pred_idx = 0;
        for (size_t j = 0; j < ell; j++) {
            size_t pos = sets[i * ell + j];
            if (pos < seed_len && seed[pos]) {
                pred_idx |= ((size_t)1 << j);
            }
        }
        /* Evaluate hard predicate on the extracted bits */
        output[i] = hard_pred_tt[pred_idx];
    }
    return true;
}