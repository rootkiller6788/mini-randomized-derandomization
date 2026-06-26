/**
 * @file nw_core.h
 * @brief Core definitions for the Nisan-Wigderson Pseudorandom Generator
 *
 * This module implements the Nisan-Wigderson (1994) construction that
 * converts a hard function into a pseudorandom generator using
 * combinatorial designs. The key insight: hardness implies randomness.
 *
 * Reference:
 *   Nisan & Wigderson, "Hardness vs Randomness", J. Comput. Syst. Sci. 49(2), 1994
 *   Arora & Barak, "Computational Complexity: A Modern Approach", Ch 20
 *
 * Notation:
 *   A (k, m, l)-design is a family S_1,...,S_m ⊆ [k] with |S_i| = l
 *   and |S_i ∩ S_j| ≤ log m for all i ≠ j.
 *   Given hard f: {0,1}^l → {0,1}, define generator G: {0,1}^k → {0,1}^m
 *   by G(x)_i = f(x|_{S_i}).
 */

#ifndef NW_CORE_H
#define NW_CORE_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ──────────────────────────────────────────────
 * L1: Core Type Definitions
 * ────────────────────────────────────────────── */

/**
 * @brief A combinatorial (k, m, l)-design.
 *
 * Universe U = {0, ..., k-1}.
 * Family of m subsets, each of size l.
 * Pairwise intersection ≤ log2(m) (or a specified bound).
 *
 * L1: Definition — NW combinatorial design
 * L3: Mathematical Structure — Set system with intersection constraints
 */
typedef struct {
    size_t k;             /**< Universe size (seed length of PRG) */
    size_t m;             /**< Number of sets (output length of PRG) */
    size_t l;             /**< Size of each set */
    size_t intersect_bound; /**< Maximum allowed pairwise intersection */
    size_t **sets;        /**< m arrays: sets[i][j] ∈ {0,...,k-1} */
} NWDesign;

/**
 * @brief A boolean function f: {0,1}^n → {0,1} with hardness guarantee.
 *
 * The hardness parameter is a lower bound on the circuit size required
 * to compute f within error (1/2 - hardness). If hardness = S, then
 * any circuit of size < S errs on at least (1/2 - hardness) fraction of inputs.
 *
 * L1: Definition — Hard function
 * L2: Core Concept — Average-case hardness
 */
typedef struct {
    size_t n;             /**< Input length in bits */
    bool (*evaluate)(const bool *input, size_t n); /**< Function pointer */
    double hardness;      /**< Circuit size lower bound */
    double epsilon;       /**< Advantage over random guessing (ε-hardness) */
    char *description;    /**< Description string */
} HardFunction;

/**
 * @brief Forward declaration of circuit type (defined in nw_circuits.h).
 *
 * L1: Definition — Boolean circuit
 * L3: Mathematical Structure — DAG representation
 */
struct Circuit;
typedef struct Circuit Circuit;

/**
 * @brief Pseudorandom generator G: {0,1}^k → {0,1}^m.
 *
 * L1: Definition — Pseudorandom Generator (PRG)
 * L2: Core Concept — Stretching (k < m)
 */
typedef struct {
    size_t seed_len;      /**< k: seed length */
    size_t output_len;    /**< m: output length */
    double stretch;       /**< m/k ratio */
    NWDesign *design;     /**< Underlying combinatorial design */
    HardFunction *hard_fn;  /**< Hard function driving the PRG */
} NWPRG;

/**
 * @brief Result of a circuit lower bound computation.
 *
 * L4: Fundamental Law — Circuit lower bounds
 */
typedef struct {
    size_t n;             /**< Number of inputs */
    size_t proven_lower_bound; /**< Proven minimum circuit size */
    double hardness;      /**< Resulting hardness parameter */
    bool is_optimal;      /**< Whether bound is known to be tight */
} CircuitBound;

/**
 * @brief Truth table for a boolean function on n variables.
 *
 * Size = 2^n bits.
 *
 * L3: Mathematical Structure — Truth table representation
 */
typedef struct {
    size_t n;             /**< Number of variables */
    size_t size;          /**< 2^n */
    bool *table;          /**< Function values in lexicographic order */
} TruthTable;

/* ──────────────────────────────────────────────
 * L2: Core Concept — Pseudorandomness
 * ────────────────────────────────────────────── */

/**
 * @brief Distinguisher: a circuit trying to tell PRG output from truly random.
 *
 * L2: Core Concept — Distinguisher / Next-bit predictor
 * L1: Definition — Computational indistinguishability
 */
typedef struct {
    size_t input_len;     /**< m: length of string to distinguish */
    Circuit *circuit;     /**< Circuit implementing the distinguisher */
    double advantage;     /**< |Pr[D(G(U_k))=1] - Pr[D(U_m)=1]| */
} Distinguisher;

/**
 * @brief Next-bit predictor: given first i bits, predicts bit i+1.
 *
 * L2: Core Concept — Equivalence of next-bit unpredictability and PRG
 * Theorem (Blum-Micali, Yao): A generator passes all polynomial-time
 * statistical tests iff it is next-bit unpredictable.
 */
typedef struct {
    size_t n;             /**< Total output length */
    size_t position;      /**< Position being predicted (0-indexed) */
    Circuit *predictor;
    double advantage;     /**< Success probability - 1/2 */
} NextBitPredictor;

/* ──────────────────────────────────────────────
 * API: Core NW Design Operations
 * ────────────────────────────────────────────── */

/**
 * @brief Create a new combinatorial design.
 *
 * @param k  Universe size (0..k-1)
 * @param m  Number of sets
 * @param l  Size of each set
 * @param intersect_bound  Maximum allowed pairwise intersection
 * @return Newly allocated NWDesign
 */
NWDesign *nw_design_create(size_t k, size_t m, size_t l, size_t intersect_bound);

/**
 * @brief Free a combinatorial design.
 * @param d  Design to free
 */
void nw_design_free(NWDesign *d);

/**
 * @brief Verify that a design satisfies the small intersection property.
 *
 * Checks ∀ i≠j: |S_i ∩ S_j| ≤ d->intersect_bound.
 * Complexity: O(m² · l)
 *
 * @param d  Design to verify
 * @return true if the intersection property holds
 */
bool nw_design_verify(const NWDesign *d);

/**
 * @brief Compute the actual maximum intersection size in a design.
 *
 * @param d  Design to analyze
 * @return max_{i≠j} |S_i ∩ S_j|
 */
size_t nw_design_max_intersection(const NWDesign *d);

/**
 * @brief Print a design in human-readable form.
 * @param d  Design to print
 * @param fp Output stream
 */
void nw_design_print(const NWDesign *d, FILE *fp);

/* ──────────────────────────────────────────────
 * API: PRG Construction and Evaluation
 * ────────────────────────────────────────────── */

/**
 * @brief Construct an NW pseudorandom generator.
 *
 * Given a (k, m, l)-design D and a hard function f: {0,1}^l → {0,1},
 * construct G: {0,1}^k → {0,1}^m where G(x)_i = f(x|_{S_i}).
 *
 * Theorem (Nisan-Wigderson 1994): If f is (S, ε)-hard and D is a
 * (k, m, l)-design with intersection ≤ log m, then G is a PRG that
 * fools circuits of size S' = S - poly(m).
 *
 * L5: Algorithm — NW PRG construction
 * L4: Fundamental Law — Nisan-Wigderson Theorem
 *
 * @param design  (k, m, l)-design
 * @param hf      Hard function f: {0,1}^l → {0,1}
 * @return Newly allocated NWPRG
 */
NWPRG *nw_prg_construct(NWDesign *design, HardFunction *hf);

/**
 * @brief Evaluate the NW PRG on a given seed.
 *
 * G(x)_i = f(x restricted to S_i) for i = 0, ..., m-1.
 *
 * Complexity: O(m · (time to evaluate f))
 *
 * @param prg    The PRG
 * @param seed   k-bit seed as bool array
 * @param output Pre-allocated output buffer of size m
 */
void nw_prg_evaluate(const NWPRG *prg, const bool *seed, bool *output);

/**
 * @brief Evaluate the NW PRG on a specific output bit.
 *
 * @param prg   The PRG
 * @param seed  k-bit seed
 * @param index Which output bit (0..m-1)
 * @return The bit value
 */
bool nw_prg_evaluate_bit(const NWPRG *prg, const bool *seed, size_t index);

/**
 * @brief Free a PRG structure.
 * @param prg  PRG to free
 */
void nw_prg_free(NWPRG *prg);

/* ──────────────────────────────────────────────
 * API: Hardness and Amplification
 * ────────────────────────────────────────────── */

/**
 * @brief Create a hard function wrapper.
 *
 * @param n          Input length
 * @param evaluate   Function pointer
 * @param hardness   Circuit-size lower bound
 * @param epsilon    Advantage over 1/2
 * @param desc       Description
 * @return Newly allocated HardFunction
 */
HardFunction *nw_hard_function_create(
    size_t n,
    bool (*evaluate)(const bool *, size_t),
    double hardness,
    double epsilon,
    const char *desc
);

/**
 * @brief Free a hard function.
 * @param hf  Hard function to free
 */
void nw_hard_function_free(HardFunction *hf);

/**
 * @brief Compute Yao's XOR Lemma amplification.
 *
 * Given f that is (1/2 - ε)-hard for circuits of size S,
 * f^{⊕k}(x_1,...,x_k) = ⊕_{i=1}^k f(x_i) is (1/2 - ε^k - δ)-hard
 * for circuits of size Ω(δ²·ε²·S).
 *
 * Theorem (Yao 1982, Goldreich-Levin 1989, Impagliazzo 1995):
 * Direct product + XOR amplifies hardness.
 *
 * L4: Fundamental Law — Yao's XOR Lemma
 * L5: Algorithm — Hardness amplification computation
 *
 * @param epsilon  Base advantage parameter
 * @param k        Number of XOR copies
 * @param delta    Slack parameter
 * @return Amplified hardness advantage
 */
double nw_yao_xor_amplify(double epsilon, size_t k, double delta);

/**
 * @brief Compute the amplified circuit size lower bound.
 *
 * @param base_size  Original circuit lower bound
 * @param epsilon    Original advantage
 * @param k          Number of XOR copies
 * @param delta      Slack
 * @return Amplified circuit lower bound
 */
double nw_amplified_circuit_size(double base_size, double epsilon, size_t k, double delta);

/* ──────────────────────────────────────────────
 * API: Seed Length Computation
 * ────────────────────────────────────────────── */

/**
 * @brief Compute the required seed length for an NW PRG.
 *
 * Given output length m and a hard function with hardness S,
 * determine the minimal seed length k such that a (k, m, l)-design
 * exists with l = O(log S).
 *
 * @param m         Desired output length
 * @param hardness  Circuit size hardness
 * @return Required seed length k
 */
size_t nw_seed_length(size_t m, double hardness);

/**
 * @brief Compute the output length achievable from a given seed.
 *
 * @param k         Seed length
 * @param hardness  Hardness parameter
 * @return Maximum output length m
 */
size_t nw_output_length(size_t k, double hardness);

/* ──────────────────────────────────────────────
 * API: Theorem Verification (Demonstration)
 * ────────────────────────────────────────────── */

/**
 * @brief Verify the NW theorem on a specific instance.
 *
 * Checks: If a next-bit predictor exists for the PRG with advantage δ,
 * then the hard function can be computed with related advantage.
 *
 * This is a constructive implementation of the hybrid argument
 * underlying the NW proof.
 *
 * L4: Fundamental Law — NW Theorem verification via hybrid argument
 *
 * @param prg         The PRG to verify
 * @param predictor   A candidate next-bit predictor
 * @param confidence  Statistical confidence level
 * @return true if the theorem's prediction holds
 */
bool nw_theorem_verify(
    const NWPRG *prg,
    const NextBitPredictor *predictor,
    double confidence
);

/**
 * @brief Run the hybrid argument: given a distinguisher for G,
 * construct a predictor for f.
 *
 * This is the core reduction in the NW proof.
 *
 * L5: Algorithm — Hybrid argument for NW PRG
 *
 * @param prg           The PRG
 * @param distinguisher A distinguisher for G(U_k) vs U_m
 * @return A predictor for f, or NULL
 */
NextBitPredictor *nw_hybrid_argument(
    const NWPRG *prg,
    const Distinguisher *distinguisher
);

/* ──────────────────────────────────────────────
 * API: Utility Operations
 * ────────────────────────────────────────────── */

/**
 * @brief Convert a uint64_t to a seed (bool array).
 * @param n      Number of bits
 * @param value  Integer value
 * @param seed   Output bool array (pre-allocated)
 */
void nw_int_to_seed(size_t n, uint64_t value, bool *seed);

/**
 * @brief Convert a bool array seed to uint64_t (for small n ≤ 64).
 * @param n     Number of bits
 * @param seed  Bool array
 * @return Integer representation
 */
uint64_t nw_seed_to_int(size_t n, const bool *seed);

/**
 * @brief Extract a subset of a seed defined by indices.
 *
 * Given seed x ∈ {0,1}^k and set S ⊆ [k],
 * return y ∈ {0,1}^{|S|} where y_j = x_{S_j}.
 *
 * This is the "restriction" operation in the NW construction.
 *
 * @param k         Original seed length
 * @param seed      Original seed
 * @param subset    Indices to extract
 * @param subset_sz Number of indices
 * @param result    Output buffer (pre-allocated, size subset_sz)
 */
void nw_seed_restrict(
    size_t k, const bool *seed,
    const size_t *subset, size_t subset_sz,
    bool *result
);

/**
 * @brief XOR two boolean arrays of equal length.
 * @param a   First array
 * @param b   Second array
 * @param out Output array
 * @param len Length
 */
void nw_bool_xor(const bool *a, const bool *b, bool *out, size_t len);

/**
 * @brief Compute Hamming weight of a boolean array.
 * @param arr  Array
 * @param len  Length
 * @return Number of true values
 */
size_t nw_hamming_weight(const bool *arr, size_t len);

/**
 * @brief Compute statistical distance between two boolean distributions
 * given as truth tables over n bits.
 *
 * @param dist1  First distribution (array of size 2^n of probabilities)
 * @param dist2  Second distribution
 * @param n      Number of variables
 * @return Total variation distance (in [0, 1])
 */
double nw_statistical_distance(const double *dist1, const double *dist2, size_t n);

/**
 * @brief Generate a truly random n-bit string (simulated).
 * @param buf  Output buffer
 * @param n    Number of bits
 */
void nw_random_bits(bool *buf, size_t n);

/**
 * @brief Compute log2 with safety for n=0.
 * @param n  Value
 * @return ceil(log2(n)), with log2(0)=0
 */
size_t nw_ceil_log2(size_t n);

#endif /* NW_CORE_H */
