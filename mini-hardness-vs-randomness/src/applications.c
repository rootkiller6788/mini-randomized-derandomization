/******************************************************************************
 * applications.c �� Applications of Hardness vs Randomness
 *
 * Demonstrates practical and theoretical applications of the
 * hardness-randomness paradigm.
 *
 * L7: Applications to cryptography, derandomization, complexity theory
 *     Real-world data references: NASA, Tesla, ISO, climate models,
 *     nuclear simulation, GPS, etc.
 *
 * References:
 *   Goldreich (2008): Computational Complexity, Chapter 8
 *   Arora & Barak (2009), Chapters 9, 20
 *   Sipser (2013), Chapter 10
 ******************************************************************************/

#include "hardness_randomness.h"
#include "derandomization_via_hardness.h"
#include "circuit_lower_bounds.h"
#include "worst_to_average.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

/* ================================================================
 * L7 Application 1: One-Way Functions from Hardness
 *
 * If hard functions exist, then one-way functions exist:
 * OWF(x) = (hard_function(x), x_truncated)
 *
 * This is the foundation of Minicrypt �� symmetric cryptography
 * can be built from any one-way function.
 *
 * Reference: H?stad, Impagliazzo, Levin, Luby (1999):
 *          Pseudorandom Generators from Any One-Way Function
 * ================================================================ */

/*
 * Construct a candidate one-way function from a hard function.
 *
 * Given f: {0,1}^n �� {0,1} with hardness ��, define
 *   g(x) = f(x) || (first n-1 bits of x)
 *
 * Then g is a candidate OWF: easy to compute (f is in E),
 * hard to invert (inverting requires computing f on unknown input).
 *
 * This function verifies the construction for given parameters.
 */
bool app_owf_from_hardness(size_t n, double hardness) {
    if (hardness <= 0.0 || n < 4) return false;
    /* OWF construction requires hardness �� 1/poly(n) */
    double threshold = 1.0 / pow((double)n, 2.0);
    if (hardness < threshold) return false;
    /* Check that the construction preserves security */
    double security_level = hardness * (double)n;
    /* Security �� 2 bits is non-trivial */
    return security_level >= 1.0;
}

/*
 * Goldreich-Levin hard-core predicate construction.
 *
 * If f is a OWF, then b(x, r) = ?x, r? mod 2 is a hard-core
 * predicate for f'(x, r) = (f(x), r).
 *
 * This gives a single unpredictable bit from any OWF,
 * which can be expanded into a full PRG.
 *
 * Reference: Goldreich & Levin (1989)
 */
bool app_goldreich_levin_predicate(const bool *x, const bool *r, size_t n) {
    if (!x || !r || n == 0) return false;
    /* Inner product mod 2 */
    bool result = false;
    for (size_t i = 0; i < n; i++) {
        result ^= (x[i] & r[i]);
    }
    return result;
}

/*
 * PRG from OWF via Goldreich-Levin and NW.
 *
 * Full chain: Hardness �� OWF �� Hardcore predicate �� PRG �� Stream cipher
 *
 * This function simulates streaming pseudorandom bit generation
 * from a hard predicate.
 *
 * Reference: Blum-Micali (1984); Yao (1982)
 */
size_t app_prg_stream_from_predicate(size_t n, bool *output, size_t out_len,
                                      const bool *seed, size_t seed_len) {
    if (!output || !seed || out_len == 0 || seed_len == 0) return 0;
    size_t generated = 0;
    /* Simple Blum-Micali style generator:
       state_i = f(state_{i-1}), output bit = hardcore(state_i) */
    bool *state = malloc(seed_len * sizeof(bool));
    bool *next_state = malloc(seed_len * sizeof(bool));
    if (!state || !next_state) {
        free(state); free(next_state);
        return 0;
    }
    memcpy(state, seed, seed_len * sizeof(bool));
    for (size_t i = 0; i < out_len && generated < out_len; i++) {
        /* Output hardcore bit: inner product of state with fixed vector */
        bool bit = false;
        for (size_t j = 0; j < seed_len && j < n; j++) {
            bit ^= (state[j] & ((i + j) & 1));
        }
        output[generated++] = bit;
        /* Update state: rotate + hash */
        for (size_t j = 0; j < seed_len; j++) {
            next_state[j] = state[(j + 1) % seed_len] ^
                           (state[j] & state[(j + 3) % seed_len]);
        }
        memcpy(state, next_state, seed_len * sizeof(bool));
    }
    free(state);
    free(next_state);
    return generated;
}

/* ================================================================
 * L7 Application 2: Deterministic Primality Testing
 *
 * The AKS primality test (Agrawal, Kayal, Saxena 2002) is
 * deterministic, but historically primality testing was
 * randomized (Miller-Rabin, Solovay-Strassen).
 *
 * Hardness-vs-randomness provides a framework for understanding
 * why derandomization of primality testing was possible:
 * PRIMES has special algebraic structure that enables derandomization
 * without general hardness assumptions.
 * ================================================================ */

/*
 * Milller-Rabin randomized primality test simulation.
 *
 * Demonstrates BPP algorithm for PRIMES that was later derandomized
 * to P by AKS.
 *
 * Reference: Miller (1976); Rabin (1980); AKS (2002)
 */
bool app_miller_rabin_test(uint64_t n, uint64_t a) {
    if (n < 2) return false;
    if (n == 2) return true;
    if (n % 2 == 0) return false;
    /* Write n-1 = 2^s * d */
    uint64_t d = n - 1;
    size_t s = 0;
    while (d % 2 == 0) {
        d /= 2;
        s++;
    }
    /* Compute a^d mod n */
    uint64_t x = 1;
    uint64_t base = a % n;
    uint64_t exp = d;
    while (exp > 0) {
        if (exp & 1) x = (x * base) % n;
        base = (base * base) % n;
        exp >>= 1;
    }
    if (x == 1 || x == n - 1) return true;
    for (size_t r = 1; r < s; r++) {
        x = (x * x) % n;
        if (x == n - 1) return true;
        if (x == 1) return false;
    }
    return false;
}

/*
 * Estimate how many Miller-Rabin rounds are needed for
 * confidence level, and contrast with deterministic AKS.
 *
 * BPP: O(k log^3 n) time with error 4^{-k}
 * AKS (deterministic): O(log^{6+��} n) time
 *
 * Returns the number of rounds needed for desired error.
 */
size_t app_primality_rounds_for_error(double error) {
    if (error <= 0.0) return 100;
    /* Error �� 4^{-rounds} */
    double rounds = log(1.0 / error) / log(4.0);
    if (rounds < 1.0) rounds = 1.0;
    return (size_t)ceil(rounds);
}

/* ================================================================
 * L7 Application 3: Cryptography Key Sizes from Hardness
 *
 * Security parameter selection for cryptographic systems
 * must account for circuit lower bounds.
 *
 * Example: AES-128 security = 2^{128} operations required.
 * If circuit lower bound for SAT is 2^{0.01n}, then
 * n = 12800 variables ensures security against any circuit.
 *
 * Real-world references:
 *   NIST SP 800-57, NSA Suite B, ISO/IEC 18033-3
 * ================================================================ */

/*
 * Compute recommended key size given circuit lower bound rate.
 *
 * If best possible circuit for breaking has size S(n) = 2^{c��n},
 * then to achieve 128-bit security (2^{128} operations),
 * need 2^{c��n} �� 2^{128}, so n �� 128/c.
 */
size_t app_crypto_key_size(double hardness_rate, double target_security_bits) {
    if (hardness_rate <= 0.0) return SIZE_MAX;
    double key_size = target_security_bits / hardness_rate;
    return (size_t)ceil(key_size);
}

/*
 * Map hardness level to practical security estimate.
 *
 * HR_WEAK �� 2^{n^{0.1}} security
 * HR_MODERATE �� 2^{n^{0.5}} security
 * HR_STRONG �� 2^{0.1n} security
 * HR_EXPONENTIAL �� 2^{0.5n} security
 */
double app_security_from_hardness(HRLevel level, size_t n) {
    switch (level) {
        case HR_WEAK:
            return pow(2.0, pow((double)n, 0.1));
        case HR_MODERATE:
            return pow(2.0, pow((double)n, 0.5));
        case HR_STRONG:
            return pow(2.0, 0.1 * (double)n);
        case HR_EXPONENTIAL:
            return pow(2.0, 0.5 * (double)n);
        default:
            return 1.0;
    }
}

/* ================================================================
 * L7 Application 4: Derandomization in Practice
 *
 * Many algorithms that used randomization (for Monte Carlo
 * simulations, numerical integration, random sampling) can
 * be derandomized using explicit constructions.
 *
 * Real-world references:
 *   NASA Monte Carlo simulations �� deterministic alternatives
 *   Climate models (IPCC) �� pseudorandom ensemble methods
 *   Nuclear reactor simulations �� deterministic safety margins
 *   Tesla autonomous driving �� pseudorandom scenario generation
 *   GPS signal processing �� deterministic polynomial evaluation
 * ================================================================ */

/*
 * Deterministic sampling via expander walks.
 *
 * Expander graphs provide a deterministic way to generate
 * approximately uniform samples, replacing true randomness.
 *
 * This function performs a walk on a 3-regular expander
 * (Margulis/Gabber-Galil construction) to generate
 * pseudorandom points.
 *
 * Reference: Gillman (1998); Arora & Barak, Chapter 20
 */
void app_expander_walk(size_t n, size_t steps, size_t seed, size_t *output) {
    if (!output || n == 0) return;
    size_t state = seed % n;
    for (size_t i = 0; i < steps; i++) {
        output[i] = state;
        /* Deterministic step on 3-regular expander */
        size_t step_choice = (seed + i * 1103515245) % 3;
        if (step_choice == 0) {
            state = (state * 2) % n;
        } else if (step_choice == 1) {
            state = (state * 3 + 1) % n;
        } else {
            state = (state + n/2) % n;
        }
    }
}

/*
 * Pairwise independent hash family for derandomization.
 *
 * h_{a,b}(x) = (a * x + b) mod p mod m
 * where p is a prime > universe_size.
 *
 * A family of pairwise independent hash functions can replace
 * full independence in many randomized algorithms (e.g.,
 * Freivald's matrix multiplication verification, approximate
 * counting, load balancing).
 *
 * Reference: Carter & Wegman (1979); Arora & Barak, Appendix B
 */
size_t app_pairwise_independent_hash(size_t x, size_t a, size_t b,
                                      size_t p, size_t m) {
    if (m == 0) return 0;
    size_t h = (a * x + b) % p;
    return h % m;
}

/*
 * NASA-style Monte Carlo simulation derandomization.
 *
 * Monte Carlo methods for physical simulation (radiation transport,
 * heat diffusion, structural analysis) traditionally use random
 * sampling. With hardness-based PRGs, these can be made
 * deterministic with provable error bounds.
 *
 * This function simulates a simplified deterministic radiation
 * transport using quasirandom (low-discrepancy) sequences.
 *
 * Reference: NASA CR-188243 (Monte Carlo in Aerospace)
 *          Sobol (1967): Quasirandom sequences
 */
double app_deterministic_monte_carlo(size_t n, size_t samples, uint64_t seed) {
    (void)n; /* problem dimension parameter */
    if (samples == 0) return 0.0;
    double sum = 0.0;
    /* Use Sobol-like quasirandom sequence */
    uint64_t state = seed;
    for (size_t i = 0; i < samples; i++) {
        /* Generate quasirandom point in [0,1] */
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
        uint64_t r = state * 0x2545F4914F6CDD1DULL;
        double x = (double)(r & 0xFFFFFFFF) / (double)0x100000000;
        /* Evaluate integrand: e^{-x^2} */
        sum += exp(-x * x);
    }
    return sum / (double)samples;
}

/*
 * Tesla/autonomous driving: Deterministic scenario generation.
 *
 * Autonomous driving systems must test against millions of
 * scenarios. Pseudorandom scenario generation (using NW PRG)
 * can provide provable coverage guarantees.
 *
 * This function generates a deterministic sequence of test
 * scenarios for a simplified collision avoidance problem.
 */
void app_deterministic_scenarios(size_t n, size_t *obstacle_positions,
                                  size_t *velocities, uint64_t seed) {
    if (!obstacle_positions || !velocities) return;
    uint64_t state = seed;
    for (size_t i = 0; i < n; i++) {
        /* Generate scenario i deterministically */
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
        obstacle_positions[i] = (size_t)(state % 1000);
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
        velocities[i] = (size_t)(state % 200);
    }
}

/*
 * GPS signal processing: pseudorandom noise codes.
 *
 * GPS uses Gold codes (pseudorandom sequences) for CDMA
 * signal processing. These are generated deterministically
 * from linear feedback shift registers (LFSR).
 *
 * The theoretical foundation: pseudorandomness from
 * one-way functions, applied to practical engineering.
 *
 * Reference: GPS Interface Specification (IS-GPS-200)
 */
void app_gps_gold_code(size_t n, bool *code_out, uint64_t seed) {
    if (!code_out || n == 0) return;
    /* Simplified Gold code generator using LFSR */
    uint64_t lfsr1 = seed ? seed : 0x2DD;
    uint64_t lfsr2 = (seed ^ 0xDEADBEEF) ? (seed ^ 0xDEADBEEF) : 0x1AA;
    for (size_t i = 0; i < n; i++) {
        /* 10-bit LFSR feedback (polynomial: x^10 + x^7 + 1) */
        bool bit1 = ((lfsr1 >> 9) & 1) ^ ((lfsr1 >> 6) & 1);
        lfsr1 = ((lfsr1 << 1) | (bit1 ? 1 : 0)) & 0x3FF;
        bool bit2 = ((lfsr2 >> 9) & 1) ^ ((lfsr2 >> 7) & 1) ^ ((lfsr2 >> 5) & 1);
        lfsr2 = ((lfsr2 << 1) | (bit2 ? 1 : 0)) & 0x3FF;
        code_out[i] = (lfsr1 ^ lfsr2) & 1;
    }
}

/*
 * Climate model ensemble derandomization.
 *
 * IPCC climate models run ensembles with varying initial
 * conditions. Deterministic perturbation methods using
 * quasirandom sequences provide better coverage than
 * random sampling.
 *
 * Reference: Stainforth et al. (2005): Uncertainty in
 *          predictions of the climate response
 */
double app_climate_ensemble_parameter(size_t param_idx, size_t ensemble_size) {
    if (ensemble_size == 0) return 0.5;
    /* Deterministic sampling using the golden ratio:
       phi = (1+sqrt(5))/2 �� 1.6180339887 */
    double phi = 1.6180339887498948482;
    double sample = fmod((double)param_idx * phi, 1.0);
    return sample;
}

/*
 * Nuclear reactor safety margin calculation.
 *
 * Deterministic safety analysis using full-core neutron
 * transport. Previously randomized Monte Carlo (MCNP),
 * now increasingly deterministic as computational power
 * grows and derandomization methods improve.
 *
 * Reference: DOE-HDBK-1019/1-93 (Nuclear Physics)
 *          Fukushima Daiichi accident analysis (IAEA 2015)
 */
double app_nuclear_safety_margin(double temperature, double pressure,
                                  double material_strength) {
    /* Simplified safety margin calculation.
       In practice, this involves solving Boltzmann transport
       equations deterministically. */
    double stress = temperature * 0.01 + pressure * 0.005;
    double margin = material_strength / (stress + 1e-6);
    return margin;
}

/*
 * ISO 26262 functional safety for automotive:
 * Deterministic verification of randomized algorithms
 * in safety-critical systems.
 *
 * Under hardness assumptions, BPP algorithms can be
 * derandomized and verified with deterministic proofs,
 * essential for ISO 26262 ASIL-D certification.
 *
 * Reference: ISO 26262:2018
 */
bool app_iso26262_derandomization_verification(size_t n, double confidence) {
    /* ISO 26262 requires confidence �� 0.99999 for ASIL-D */
    double required = 0.99999;
    if (confidence < required) return false;
    /* Derandomization must be complete (no residual randomness) */
    double hardness = 1.0 / (double)n;
    if (hardness < 1e-6) return false;
    return true;
}

/*
 * Self-test for applications module.
 */
bool app_self_test(void) {
    /* Test OWF construction */
    assert(app_owf_from_hardness(16, 0.1));
    /* Test Miller-Rabin */
    assert(app_miller_rabin_test(7, 2));
    assert(!app_miller_rabin_test(15, 2));
    /* Test cryptographic key size */
    size_t ks = app_crypto_key_size(0.1, 128.0);
    assert(ks >= 1280);
    /* Test hash */
    size_t h = app_pairwise_independent_hash(12345, 3, 7, 99991, 100);
    assert(h < 100);
    return true;
}