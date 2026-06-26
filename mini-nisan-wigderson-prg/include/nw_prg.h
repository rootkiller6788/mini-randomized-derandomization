/**
 * @file nw_prg.h
 * @brief Nisan-Wigderson Pseudorandom Generator — Full Implementation.
 *
 * This module implements the core Nisan-Wigderson (1994) construction
 * and its derandomization consequences:
 *
 *   Hardness of f  +  Combinatorial design  ⟹  Pseudorandom generator
 *   PRG with stretch  ⟹  Derandomization of BPP
 *
 * The key insight: if there exists a function in E that requires
 * exponential-size circuits, then randomness can be eliminated
 * from polynomial-time algorithms with only a sub-exponential
 * slowdown.
 *
 * Reference:
 *   Nisan & Wigderson, "Hardness vs Randomness", JCSS 49(2), 1994
 *   Impagliazzo & Wigderson, "P = BPP if E requires exponential circuits", STOC 1997
 *   Arora & Barak, "Computational Complexity", Ch 20
 *   Goldreich, "Computational Complexity: A Conceptual Perspective", Ch 8
 */

#ifndef NW_PRG_H
#define NW_PRG_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include "nw_core.h"
#include "nw_designs.h"
#include "nw_hardness.h"
#include "nw_circuits.h"

/* ──────────────────────────────────────────────
 * L2: PRG Construction Modes
 * ────────────────────────────────────────────── */

/**
 * @brief PRG construction strategy.
 *
 * L2: Core Concept — Construction method taxonomy
 */
typedef enum {
    NW_PRG_STANDARD,        /**< Standard NW: hard fn + design */
    NW_PRG_YAO_XOR,        /**< NW with XOR-amplified hardness */
    NW_PRG_IW97,           /**< IW97: full derandomization pipeline */
    NW_PRG_OPTIMIZED       /**< Optimized parameters (tight bounds) */
} NWPRGMode;

/* ──────────────────────────────────────────────
 * L3: PRG Security Analysis
 * ────────────────────────────────────────────── */

/**
 * @brief Security analysis result for an NW PRG instance.
 *
 * L4: Fundamental Law — PRG security analysis
 */
typedef struct {
    double distinguishing_advantage;  /**< ε: max advantage of distinguisher */
    double circuit_size_bound;        /**< S: circuits of size < S are fooled */
    double seed_stretch;              /**< output/seed ratio */
    bool is_secure;                   /**< true if provably secure */
    size_t security_bits;             /**< log2(circuit_size_bound) */
} PRGSecurity;

/**
 * @brief Result of running a statistical test on PRG output.
 *
 * L5: Algorithm — Statistical testing
 */
typedef struct {
    size_t num_samples;         /**< Number of independent samples tested */
    size_t num_tests;           /**< Number of statistical tests applied */
    double *p_values;           /**< p-value for each test */
    bool passes_all;            /**< true if no test rejects at α = 0.01 */
    double min_p_value;         /**< Minimum p-value observed */
} StatisticalTestResult;

/* ──────────────────────────────────────────────
 * L4: Core NW PRG Construction Algorithm
 * ────────────────────────────────────────────── */

/**
 * @brief Construct the optimal NW PRG for given parameters.
 *
 * This is the main entry point. It:
 *   1. Determines if a design exists with the required parameters
 *   2. Constructs the design (Reed-Solomon or random)
 *   3. Optionally amplifies hardness via Yao's XOR Lemma
 *   4. Sets up the PRG structure
 *
 * Theorem (Nisan-Wigderson 1994):
 * Let S: N → N. Suppose there exists f ∈ EXP such that for all
 * sufficiently large n, Circuits(S(n)) cannot compute f_n.
 * Then for any polynomial m(n), there exists a PRG
 * G: {0,1}^{O(log n)} → {0,1}^{m(n)} that fools circuits of
 * size n. In particular, BPP ⊂ SUBEXP.
 *
 * L4: Fundamental Law — Nisan-Wigderson Theorem (full statement)
 * L5: Algorithm — Optimal PRG construction
 *
 * @param output_len  m: desired output length
 * @param hardness    S: circuit-size hardness of underlying function
 * @param mode        Construction strategy
 * @return Constructed PRG, or NULL if infeasible
 */
NWPRG *nw_prg_construct_optimal(
    size_t output_len, double hardness, NWPRGMode mode
);

/**
 * @brief Compute the stretch achievable for given hardness.
 *
 * Given hardness S, the maximum output length m for a PRG with seed
 * length k is approximately:
 *   m = 2^{k^{1/3}} for sub-exponential hardness
 *   m = 2^{Ω(k)} for exponential hardness (IW97)
 *
 * L4: Fundamental Law — PRG stretch vs hardness tradeoff
 *
 * @param seed_len  k
 * @param hardness  S
 * @return Maximum output length m
 */
size_t nw_prg_max_stretch(size_t seed_len, double hardness);

/**
 * @brief Compute the minimal seed length for given output length.
 *
 * Inverts the stretch computation.
 *
 * @param output_len  m
 * @param hardness    S
 * @return Minimal seed length k
 */
size_t nw_prg_min_seed(size_t output_len, double hardness);

/* ──────────────────────────────────────────────
 * L5: PRG Evaluation and Round Operations
 * ────────────────────────────────────────────── */

/**
 * @brief Evaluate the NW PRG and collect all output bits.
 *
 * For each i = 0, ..., m-1:
 *   1. Extract x|_{S_i} from seed
 *   2. Compute f(x|_{S_i})
 *   3. Store as output[i]
 *
 * Complexity: O(m · T_f(l))
 *
 * @param prg   The PRG
 * @param seed  Seed (size k)
 * @param out   Output buffer (size m), pre-allocated
 */
void nw_prg_eval_all(const NWPRG *prg, const bool *seed, bool *out);

/**
 * @brief Evaluate the PRG on a range of output bits.
 *
 * Useful for streaming/lazy evaluation.
 *
 * @param prg         The PRG
 * @param seed        Seed
 * @param start_idx   First output bit index
 * @param count       Number of output bits
 * @param out         Output buffer, pre-allocated
 */
void nw_prg_eval_range(
    const NWPRG *prg, const bool *seed,
    size_t start_idx, size_t count, bool *out
);

/**
 * @brief Evaluate PRG with 64-bit word seed representation.
 *
 * @param prg       The PRG
 * @param seed_word 64-bit seed (lower k bits used)
 * @param out       Output buffer (size m)
 */
void nw_prg_eval_word(const NWPRG *prg, uint64_t seed_word, bool *out);

/**
 * @brief Evaluate PRG with 64-bit output (for m ≤ 64).
 *
 * @param prg       The PRG
 * @param seed_word 64-bit seed
 * @return 64-bit output (lower m bits)
 */
uint64_t nw_prg_eval_word_to_word(const NWPRG *prg, uint64_t seed_word);

/* ──────────────────────────────────────────────
 * L5: Distinguisher Simulation
 * ────────────────────────────────────────────── */

/**
 * @brief Run a distinguisher circuit on PRG output vs random.
 *
 * This simulates the security analysis: given a candidate distinguisher
 * circuit D, estimate its advantage ε = |Pr[D(G(U_k))=1] - Pr[D(U_m)=1]|.
 *
 * L5: Algorithm — Distinguisher advantage estimation
 *
 * @param prg              The PRG
 * @param distinguish_circuit  Circuit D: {0,1}^m → {0,1}
 * @param num_samples      Number of Monte Carlo samples
 * @return Estimated advantage
 */
double nw_prg_distinguisher_advantage(
    const NWPRG *prg,
    const Circuit *distinguish_circuit,
    size_t num_samples
);

/**
 * @brief Construct a next-bit predictor from a distinguisher.
 *
 * This is the constructive version of Yao's theorem: if a distinguisher
 * D can tell G(U_k) from U_m with advantage ε, then there exists a
 * next-bit predictor P_i that predicts the (i+1)-th bit given the first
 * i bits with advantage ε/m.
 *
 * L4: Fundamental Law — Yao's distinguisher-to-predictor theorem
 * L5: Algorithm — Predictor construction via hybrid argument
 *
 * @param prg              The PRG
 * @param distinguisher    Circuit D
 * @param advantage        ε
 * @return Array of m next-bit predictors (one per position)
 */
NextBitPredictor *nw_distinguisher_to_predictors(
    const NWPRG *prg,
    const Circuit *distinguisher,
    double advantage
);

/* ──────────────────────────────────────────────
 * L7: BPP Derandomization via NW PRG
 * ────────────────────────────────────────────── */

/**
 * @brief Derandomize a BPP algorithm using the NW PRG.
 *
 * Given:
 *   A(x, r): a randomized algorithm for language L with error ≤ 1/3
 *   where |r| = m(n) random bits are used on inputs of length n
 *   NW PRG G: {0,1}^k → {0,1}^m
 *
 * Output: deterministic algorithm that on input x enumerates all
 * seeds s ∈ {0,1}^k, runs A(x, G(s)), and takes majority vote.
 *
 * Complexity: 2^k · (T_A(n) + T_G(n))
 *
 * L7: Application — BPP derandomization
 * L5: Algorithm — BPP simulation via NW PRG
 *
 * @param prg         NW PRG
 * @param bp_algorithm  The BPP algorithm (takes n-bit input, m-bit random)
 * @param input       Input to the algorithm
 * @param n           Input length
 * @param m           Number of random bits used
 * @param result      Output: result of majority vote
 * @return true if exactly one output had majority
 */
bool nw_bpp_derandomize(
    const NWPRG *prg,
    bool (*bp_algorithm)(const bool *input, const bool *random, size_t n, size_t m),
    const bool *input, size_t n, size_t m,
    bool *result
);

/**
 * @brief Construct a PRG suitable for derandomizing BPP algorithms
 * of input size n using sub-exponential hardness.
 *
 * L7: Application — Sub-exponential derandomization
 *
 * @param n                Input size
 * @param randomness       Number of random bits used by BPP algorithm
 * @param hardness         Hardness parameter (> 0)
 * @return Constructed PRG
 */
NWPRG *nw_prg_for_BPP(size_t n, size_t randomness, double hardness);

/* ──────────────────────────────────────────────
 * L5: Statistical Testing Suite
 * ────────────────────────────────────────────── */

/**
 * @brief Run the NIST-like statistical test suite on PRG output.
 *
 * Tests:
 *   1. Monobit (frequency) test
 *   2. Runs test
 *   3. Longest run of ones test
 *   4. Binary matrix rank test
 *   5. Discrete Fourier transform test
 *
 * This is for demonstration/testing purposes; the NW PRG's security
 * is based on computational indistinguishability, not statistical tests.
 *
 * L5: Algorithm — Statistical testing of PRG output
 *
 * @param prg          The PRG
 * @param num_streams  Number of independent streams to test
 * @param stream_len   Length of each stream in bits
 * @return Test results
 */
StatisticalTestResult *nw_prg_statistical_tests(
    const NWPRG *prg, size_t num_streams, size_t stream_len
);

/**
 * @brief Free a statistical test result.
 * @param r  Result to free
 */
void nw_statistical_test_result_free(StatisticalTestResult *r);

/* ──────────────────────────────────────────────
 * L5: Hybrid Argument Implementation
 * ────────────────────────────────────────────── */

/**
 * @brief Implement the full hybrid argument for NW PRG.
 *
 * Hybrid H_i: first i bits from PRG on random seed, remaining
 * m-i bits truly random.
 * By hybrid argument, if |Pr[D(H_i)=1] - Pr[D(H_{i-1})=1]| ≥ ε/m
 * for some i, then we can predict the i-th PRG output bit from
 * the first i-1 bits.
 *
 * L5: Algorithm — Hybrid argument for cryptographic reductions
 *
 * @param prg            The PRG
 * @param distinguisher  Circuit D
 * @param num_samples    Samples per hybrid
 * @return Index i with largest gap, or -1 if none significant
 */
int nw_hybrid_position(
    const NWPRG *prg,
    const Circuit *distinguisher,
    size_t num_samples
);

/* ──────────────────────────────────────────────
 * L7: Cryptographic Applications
 * ────────────────────────────────────────────── */

/**
 * @brief Use the NW PRG as a key derivation function.
 *
 * Given a short random seed (key), expand to a longer
 * pseudorandom keystream.
 *
 * L7: Application — Cryptography (stream cipher from PRG)
 *
 * @param prg        The PRG
 * @param seed       Short random key
 * @param keystream  Output buffer for derived bits
 * @param stream_len Number of bits needed
 */
void nw_prg_key_derivation(
    const NWPRG *prg, const bool *seed,
    bool *keystream, size_t stream_len
);

/**
 * @brief Use the NW PRG for randomized algorithm simulation
 * with rigorous error bounds.
 *
 * L7: Application — Randomized algorithm derandomization
 *
 * @param prg          The PRG
 * @param algorithm    Randomized algorithm to simulate
 * @param input        Algorithm input data
 * @param input_size   Size of input in bits
 * @param rand_bits    Number of random bits algorithm expects
 * @param error_bound  Target error probability
 * @return Deterministic output
 */
bool nw_prg_simulate_randomized(
    const NWPRG *prg,
    bool (*algorithm)(const void *input, const bool *random, size_t in_sz, size_t r_sz),
    const void *input, size_t input_size,
    size_t rand_bits, double error_bound
);

/* ──────────────────────────────────────────────
 * Utility: PRG State Management
 * ────────────────────────────────────────────── */

/**
 * @brief Compute security level (in bits) achieved by the PRG.
 *
 * @param prg  The PRG
 * @return Security level
 */
size_t nw_prg_security_level(const NWPRG *prg);

/**
 * @brief Print PRG parameters for debugging/inspection.
 * @param prg  The PRG
 * @param fp   Output stream
 */
void nw_prg_print(const NWPRG *prg, FILE *fp);

/**
 * @brief Clone a PRG (deep copy).
 * @param prg  Source PRG
 * @return Deep copy
 */
NWPRG *nw_prg_clone(const NWPRG *prg);

/**
 * @brief Default hard function: inner product mod 2.
 *
 * IP(x_1,...,x_n) = XOR_{i=1}^{n/2} (x_i AND x_{i+n/2})
 *
 * This is the Goldreich-Levin hard-core predicate, provably hard
 * under the assumption that one-way functions exist.
 *
 * L2: Core Concept — Candidate hard function for NW PRG
 *
 * @param input  n-bit input
 * @param n      Input length
 * @return Inner product parity
 */
bool nw_default_hard_function(const bool *input, size_t n);

#endif /* NW_PRG_H */
