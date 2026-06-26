/******************************************************************************
 * iw_theorem.c �� Impagliazzo-Wigderson Theorem and Five Worlds
 *
 * Implements the full Impagliazzo-Wigderson framework:
 * BPP = P if there exists a function in EXP with exponential
 * circuit complexity (E requires circuits of size 2^{��(n)}).
 *
 * Also covers Impagliazzo's "Five Worlds" of average-case complexity.
 *
 * Knowledge: L1 (P, BPP, E, EXP definitions), L2 (IW theorem concept),
 * L4 (IW theorem proof structure), L7 (cryptographic implications),
 * L8 (Impagliazzo's five worlds)
 *
 * References:
 *   Impagliazzo & Wigderson (1997): P = BPP if E requires exponential circuits
 *   Impagliazzo (1995): A personal view of average-case complexity
 *   Arora & Barak (2009), Chapter 20
 ******************************************************************************/

#include "hardness_randomness.h"
#include "derandomization_via_hardness.h"
#include "circuit_lower_bounds.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

/* ================================================================
 * L1: Impagliazzo's Five Worlds
 *
 * Impagliazzo (1995) described five possible worlds for the
 * relationship between worst-case and average-case complexity:
 *
 * Algorithmica:  P = NP (or BPP = P, no cryptography)
 * Heuristica:   NP hard in worst case, but easy on average
 * Pessiland:    NP hard on average, but no one-way functions
 * Minicrypt:    One-way functions exist (and hence symmetric crypto)
 * Cryptomania:  Public-key cryptography exists
 *
 * The hardness-vs-randomness paradigm provides evidence that
 * we live in Minicrypt or Cryptomania (hardness implies PRG,
 * PRG implies OWF, OWF implies crypto).
 * ================================================================ */

/* The five worlds enumeration */
typedef enum {
    WORLD_ALGORITHMICA = 0,
    WORLD_HEURISTICA,
    WORLD_PESSILAND,
    WORLD_MINICRYPT,
    WORLD_CRYPTOMANIA
} ImpagliazzoWorld;

/*
 * Determine which world we are in based on hardness assumptions.
 *
 * Algorithmica:  No strong hardness found �� P = NP plausible
 * Heuristica:    Worst-case hard, but average-case easy
 * Pessiland:     Average hard but no OWF (requires extra conditions)
 * Minicrypt:     Hardness implies OWF �� symmetric crypto exists
 * Cryptomania:   Strong OWF �� PKE exists (trapdoor functions)
 *
 * Returns a descriptive string for the world.
 */
const char *iw_world_name(ImpagliazzoWorld world) {
    switch (world) {
        case WORLD_ALGORITHMICA: return "Algorithmica";
        case WORLD_HEURISTICA:  return "Heuristica";
        case WORLD_PESSILAND:   return "Pessiland";
        case WORLD_MINICRYPT:   return "Minicrypt";
        case WORLD_CRYPTOMANIA: return "Cryptomania";
        default:                return "Unknown";
    }
}

/*
 * Classify the complexity world from hardness parameters.
 *
 * Input: worst-case hardness �� and average-case hardness ��.
 * The classification follows Impagliazzo's logic:
 *
 * - If �� is negligible (no hard functions) �� Algorithmica
 * - If �� > 0 but �� is small �� Heuristica
 * - If both �� and �� > 0 but no OWF construction �� Pessiland
 * - If hardness yields OWF �� Minicrypt
 * - If trapdoor permutation exists �� Cryptomania
 */
ImpagliazzoWorld iw_classify_world(double worst_hardness, double avg_hardness,
                                    bool owf_from_hardness, bool trapdoor) {
    if (worst_hardness < 0.001) {
        return WORLD_ALGORITHMICA;
    }
    if (avg_hardness < 0.001) {
        return WORLD_HEURISTICA;
    }
    if (!owf_from_hardness) {
        return WORLD_PESSILAND;
    }
    if (trapdoor) {
        return WORLD_CRYPTOMANIA;
    }
    return WORLD_MINICRYPT;
}

/*
 * IW Theorem: BPP = P under exponential hardness.
 *
 * If there exists f in E (DTIME(2^{O(n)})) with circuit complexity
 * 2^{��(n)}, then:
 *   1. BPP ? P (derandomization)
 *   2. There exists PRG that fools polynomial-size circuits
 *   3. The NW construction yields deterministic simulations
 *
 * This function checks whether the IW conditions are met
 * for a given hardness certificate and input size.
 *
 * Reference: Impagliazzo & Wigderson (1997), Theorem 1
 */
bool iw_theorem_conditions_met(const HardnessCert *cert, size_t n) {
    if (!cert || n == 0) return false;
    /* Need: hardness �� c/n for some constant c */
    double min_hardness = 0.5 / (double)n;
    if (cert->hardness < min_hardness) return false;
    /* Need: lower bound is against circuits of size �� n */
    if (cert->circuit_sz < 2) return false;
    return true;
}

/*
 * Compute the IW derandomization time.
 *
 * Under IW assumptions, BPP can be simulated in deterministic time:
 *   T_det(n) = n^{O(1)} * T_BPP(n)
 *
 * The exponent depends on the hardness:
 *   time_blowup = O(n^{c/hardness})
 *
 * Reference: Impagliazzo & Wigderson (1997)
 */
double iw_derandomization_blowup(size_t n, double hardness) {
    if (hardness <= 0.0) return HUGE_VAL;
    /* Blowup factor = n^{O(1/hardness)} */
    double exp = 2.0 / hardness;
    return pow((double)n, exp);
}

/*
 * Exponential Complexity Assumption (ECA):
 *
 * There exists f in E = DTIME(2^{O(n)}) such that
 * for every family of circuits {C_n} of size �� 2^{�š�n},
 * C_n fails to compute f_n on at least 1/2 - 2^{-�š�n} inputs.
 *
 * This is the standard assumption that implies BPP = P.
 *
 * This function tests whether ECA is plausible for given
 * parameters by checking consistency.
 */
bool iw_check_eca_assumption(size_t n, double epsilon) {
    if (n == 0 || epsilon <= 0.0 || epsilon > 1.0) return false;
    /* ECA: need �� > 0 for the assumption to be meaningful */
    if (epsilon < 0.001) return false;
    /* Circuit size bound: 2^{�š�n} must be less than 2^{2^n}
       (which is the total number of functions) */
    double circuit_bound = pow(2.0, epsilon * (double)n);
    double total_functions = pow(2.0, pow(2.0, (double)n));
    if (circuit_bound >= total_functions) return false;
    return true;
}

/*
 * Compute the circuit lower bound threshold needed for IW.
 *
 * For the IW theorem to yield BPP = P at input length n,
 * we need a function in E with circuit complexity �� S_E(n)
 * where S_E(n) = 2^{��(n)}.
 *
 * This function computes S_E(n) for given n.
 */
size_t iw_required_circuit_size(size_t n) {
    if (n < 2) return 2;
    /* S_E(n) = 2^{c��n} for c �� 0.1 (conservative) */
    return (size_t)pow(2.0, 0.1 * (double)n);
}

/*
 * IW Full Derandomization Theorem:
 *
 * Theorem (Impagliazzo-Wigderson 1997):
 * If DTIME(2^{O(n)}) contains a language requiring circuits
 * of size 2^{��(n)}, then P = BPP.
 *
 * Proof structure:
 * 1. Hardness amplification: worst-case �� average-case (via XOR lemma)
 * 2. NW construction: average-case hard function �� PRG
 * 3. PRG-based derandomization: replace random bits with
 *    pseudorandom bits, enumerate all seeds
 *
 * This function simulates the full pipeline for small n.
 */
bool iw_full_derandomization_pipeline(size_t n, HRLevel level) {
    if (n > 20) return false;
    /* Step 1: Check if hardness level provides sufficient lower bound */
    double min_lb = clb_min_lower_bound_for_level(level, n);
    size_t circuit_lb = (size_t)min_lb;
    /* Step 2: Create hardness certificate */
    HardnessCert cert = derand_cert_from_circuit_lb(n, circuit_lb);
    /* Step 3: Verify certificate */
    if (!derand_verify_certificate(&cert, n)) return false;
    /* Step 4: Check IW conditions */
    if (!iw_theorem_conditions_met(&cert, n)) return false;
    /* Step 5: Compute PRG parameters */
    size_t seed_len = hr_seed_len(n, level);
    size_t output_len = hr_output_len(n, level);
    /* Step 6: Verify PRG stretch is meaningful */
    if (output_len <= seed_len) return false;
    /* All conditions met �� derandomization is possible */
    return true;
}

/*
 * Bridge from HR to BPP derandomization.
 *
 * This function integrates the hardness_randomness module
 * with the derandomization module to provide a unified API
 * for checking whether BPP = P under given assumptions.
 */
bool iw_bpp_subset_p(size_t n, HRLevel level) {
    /* Check derandomization feasibility */
    bool feasible = hr_bpp_equals_p_check(level, n);
    if (!feasible) return false;
    /* Compute required parameters */
    double hardness = hr_hardness_from_circuit_lb(n,
        (size_t)clb_min_lower_bound_for_level(level, n));
    /* Verifiy hardness is sufficient */
    if (hardness < hr_impagliazzo_wigderson_threshold(n)) return false;
    return true;
}

/*
 * Exponential Time Hypothesis (ETH) version for IW:
 *
 * If 3SAT requires 2^{��(n)} time, then by the Cook-Levin
 * reduction, there are hard functions in E, which via IW
 * implies BPP = P.
 *
 * This computes the gap between ETH and the required
 * hardness for derandomization.
 */
double iw_eth_derandomization_gap(size_t n) {
    if (n < 2) return 0.0;
    /* ETH: 3SAT time �� 2^{c��n}, c �� 0.05 */
    double eth_hardness = 0.05;
    /* IW required hardness */
    double iw_threshold = hr_impagliazzo_wigderson_threshold(n);
    /* Gap = how much stronger than needed */
    return eth_hardness - iw_threshold;
}

/*
 * Non-black-box derandomization bound.
 *
 * The IW theorem is inherently non-black-box: it uses the
 * circuit representation of the BPP algorithm, not just
 * oracle access. This function computes bounds on the
 * derandomization efficiency in the non-black-box setting.
 *
 * Reference: Impagliazzo & Wigderson (1997), Section 3
 */
size_t iw_non_black_box_time(size_t n, double hardness) {
    if (hardness <= 0.0) return SIZE_MAX;
    /* Non-black-box: time = poly(n) * 2^{O(log n / hardness)} */
    double exponent = 3.0 * log2((double)n) / hardness;
    double time = pow(2.0, exponent);
    return (size_t)time;
}

/*
 * Average-case derandomization under NW assumptions.
 *
 * Instead of full BPP = P, we can get average-case
 * derandomization: for most inputs, the deterministic
 * simulation works correctly.
 *
 * References: Nisan & Wigderson (1994), Section 5
 */
bool iw_average_case_derandomization_check(size_t n, double hardness) {
    if (hardness <= 0.0) return false;
    /* Average-case derandomization requires slightly weaker
       assumptions: hardness �� poly(log n)/n */
    double threshold = log2((double)n + 1) / (double)n;
    return hardness >= threshold;
}

/*
 * Pseudo-determinism from hardness:
 *
 * Under hardness assumptions, BPP algorithms can be made
 * pseudo-deterministic: they produce the same output on
 * most random choices.
 *
 * Reference: Gat & Goldwasser (2011); Arora & Barak, Chapter 20
 */
double iw_pseudodeterministic_probability(size_t n, double hardness) {
    if (hardness <= 0.0) return 0.0;
    /* Probability �� 1 - 2^{-n} for strong hardness */
    double prob = 1.0 - pow(2.0, -(double)n * hardness);
    if (prob < 0.0) prob = 0.0;
    if (prob > 1.0) prob = 1.0;
    return prob;
}

/*
 * Circuit lower bound check for IW construction.
 *
 * The NW construction requires a function that is hard
 * against circuits of size n^{O(1)}. This checks whether
 * the given lower bound is sufficient.
 */
bool iw_circuit_lb_sufficient(const CircuitLB *clb) {
    if (!clb) return false;
    /* Need: lower_bound > n^2 (at minimum for non-trivial) */
    if (clb->lower_bound < (double)(clb->n * clb->n)) return false;
    /* Need: lower bound is exponential in n */
    double quality = clb_quality_metric(clb);
    /* Quality = log(lower_bound)/n. Need > 0.1 for meaningful IW. */
    return quality >= 0.1;
}

/*
 * Self-test: verify IW theorem pipeline with known examples.
 *
 * Tests the pipeline with PARITY lower bounds (AC0) as a
 * demonstration case (even though PARITY is not in E with 
 * exponential circuit lower bounds �� it's in AC0-complete
 * for weak circuit classes).
 */
bool iw_self_test(void) {
    /* Test with moderate parameters */
    size_t n = 16;
    CircuitLB lb = clb_shannon_counting(n);
    assert(lb.lower_bound > 100.0);
    assert(clb_is_strong_enough(&lb, HR_STRONG));
    HardnessCert cert = derand_cert_from_circuit_lb(n, (size_t)lb.lower_bound);
    assert(cert.hardness > 0.0);
    bool ok = iw_theorem_conditions_met(&cert, n);
    return ok;
}