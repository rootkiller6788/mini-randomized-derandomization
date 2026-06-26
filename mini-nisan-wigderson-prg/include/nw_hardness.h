/**
 * @file nw_hardness.h
 * @brief Hardness assumptions, amplification, and truth table analysis.
 *
 * The Nisan-Wigderson construction requires a function that is hard
 * on average for small circuits. This module provides:
 *   1. Hardness assumption types and validation
 *   2. Yao's XOR Lemma for hardness amplification
 *   3. Worst-case to average-case reduction (Impagliazzo 1995)
 *   4. Truth table representation and hardness testing
 *   5. Natural proofs barrier analysis
 *
 * Reference:
 *   Yao, "Theory and Applications of Trapdoor Functions", FOCS 1982
 *   Impagliazzo, "Hard-core distributions for somewhat hard problems", FOCS 1995
 *   Impagliazzo & Wigderson, "P = BPP if E requires exponential circuits", STOC 1997
 *   Razborov & Rudich, "Natural Proofs", J. Comput. Syst. Sci. 55(1), 1997
 */

#ifndef NW_HARDNESS_H
#define NW_HARDNESS_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include "nw_core.h"

/* ──────────────────────────────────────────────
 * L1: Hardness Assumption Types
 * ────────────────────────────────────────────── */

/**
 * @brief A hardness assumption: E requires exponential-size circuits.
 *
 * This is the standard cryptographic hardness assumption.
 * Variants:
 *   - E-hard: ∃f ∈ DTIME(2^{O(n)}) requiring circuits of size 2^{Ω(n)}
 *   - Subexp-hard: circuits of size 2^{n^ε} for some ε > 0
 *
 * L1: Definition — Hardness assumption
 * L2: Core Concept — Circuit lower bound as hardness
 */
typedef struct {
    size_t n;             /**< Input length */
    double circuit_size;  /**< Required circuit size lower bound */
    double hardness;      /**< Derived hardness parameter */
    bool is_exponential;  /**< true = 2^{Ω(n)}, false = 2^{n^o(1)} */
} HardnessAssumption;

/**
 * @brief Worst-case to average-case reduction state.
 *
 * Impagliazzo (1995): If f is δ-hard in the worst case for circuits
 * of size S, then there exists a "hard-core" set of density δ on
 * which f is (1/2 - ε)-hard for circuits of size poly(ε, δ) · S.
 *
 * L4: Fundamental Law — Impagliazzo's hard-core set theorem
 * L5: Algorithm — Hard-core distribution construction
 */
typedef struct {
    size_t n;             /**< Input length */
    double worst_hardness;  /**< Worst-case hardness parameter */
    double avg_hardness;    /**< Average-case hardness achieved */
    double hardcore_density; /**< Fraction of inputs in hard-core set */
    bool *hardcore_set;     /**< Indicator for hard-core inputs (size 2^n) */
} WorstToAverage;

/**
 * @brief Compression-based representation for hard functions.
 *
 * L8: Advanced Topic — Minimal circuit representation
 */
typedef struct {
    size_t n;             /**< Number of inputs */
    size_t compressed_sz; /**< Compressed representation size */
    uint8_t *data;        /**< Compressed function data */
} CompressedFunction;

/* ──────────────────────────────────────────────
 * L4/L5: Yao's XOR Lemma
 * ────────────────────────────────────────────── */

/**
 * @brief Apply Yao's XOR Lemma to a hard function.
 *
 * Given f: {0,1}^n → {0,1} that is (ε, S)-hard (i.e., no circuit
 * of size S can compute f on > 1-ε fraction of inputs), XOR-k
 * produces g: {0,1}^{k·n} → {0,1} that is (ε^k + δ, S')-hard
 * with S' = Ω(δ² · ε² · S) for any δ > 0.
 *
 * Theorem (Yao 1982, Levin 1987, Goldreich et al. 2011):
 * The XOR lemma amplifies hardness exponentially in k.
 *
 * L4: Fundamental Law — Yao's XOR Lemma (full statement)
 * L5: Algorithm — XOR amplification computation
 *
 * @param n          Input length of base function
 * @param epsilon    Base hardness parameter (advantage)
 * @param k          Number of independent copies to XOR
 * @param delta      Slack parameter (> 0)
 * @param base_size  Base circuit lower bound S
 * @param[out] new_n       Output: new input length (k·n)
 * @param[out] new_eps     Output: amplified advantage
 * @param[out] new_size    Output: amplified circuit lower bound
 */
void nw_yao_xor_lemma_full(
    size_t n, double epsilon, size_t k, double delta, double base_size,
    size_t *new_n, double *new_eps, double *new_size
);

/**
 * @brief Compute the XOR of k truth table entries.
 *
 * For a truth table tt of size 2^n, compute:
 *   g(x_1, ..., x_k) = ⊕_{i=1}^k f(x_i)
 * where each x_i is an n-bit string.
 *
 * @param tt  Base truth table
 * @param k   Number of copies
 * @return XOR-amplified truth table (size 2^{k·n}), caller must free
 */
TruthTable *nw_truth_table_xor_k(const TruthTable *tt, size_t k);

/* ──────────────────────────────────────────────
 * L4: Worst-Case to Average-Case Reduction
 * ────────────────────────────────────────────── */

/**
 * @brief Convert worst-case hardness to average-case hardness.
 *
 * Impagliazzo's hard-core set theorem: If f is δ-hard in the worst
 * case (any circuit of size S fails on at least δ fraction of inputs),
 * then there is a set H of density ≥ δ such that f is (1/2 - ε)-hard
 * on H for circuits of size ε²·δ²·S.
 *
 * This function computes the hard-core set via the boosting algorithm
 * of Klivans & Servedio (2003).
 *
 * L4: Fundamental Law — Impagliazzo's hard-core lemma
 * L5: Algorithm — Hard-core set construction
 *
 * @param tt             Truth table of the function
 * @param worst_hardness δ (fraction of hard inputs)
 * @param epsilon        Target average-case advantage
 * @return WorstToAverage structure with hard-core set
 */
WorstToAverage *nw_worst_to_average(
    const TruthTable *tt, double worst_hardness, double epsilon
);

/**
 * @brief Free a WorstToAverage structure.
 * @param wta  Structure to free
 */
void nw_worst_to_average_free(WorstToAverage *wta);

/* ──────────────────────────────────────────────
 * L5: Hardness Assumption Utilities
 * ────────────────────────────────────────────── */

/**
 * @brief Validate a hardness assumption.
 *
 * Checks:
 *   1. circuit_size > 0
 *   2. hardness ∈ (0, 1/2]
 *   3. If exponential, circuit_size ≥ 2^{n·c} for c > 0
 *
 * @param ha  Assumption to validate
 * @return true if the assumption is well-formed
 */
bool nw_hardness_assumption_valid(const HardnessAssumption *ha);

/**
 * @brief Compute the circuit size lower bound for a given n.
 *
 * Uses Lupanov's bound: almost all functions require 2^n/n gates.
 * Shannon (1949): ∃f requiring 2^n/(2n) gates.
 *
 * L4: Fundamental Law — Shannon's counting argument
 *
 * @param n     Number of inputs
 * @param gates Actual gate count of some candidate circuit
 * @return Minimum possible circuit size for a non-trivial function
 */
double nw_shannon_lower_bound(size_t n);

/**
 * @brief Compute the Lupanov upper bound: every function has
 * a circuit of size (1 + o(1)) 2^n/n.
 *
 * @param n  Number of inputs
 * @return Upper bound on circuit complexity
 */
double nw_lupanov_upper_bound(size_t n);

/* ──────────────────────────────────────────────
 * L5: Truth Table Operations
 * ────────────────────────────────────────────── */

/**
 * @brief Create a truth table from a function pointer.
 *
 * Evaluates f on all 2^n inputs in lexicographic order.
 * Complexity: O(2^n · T_f(n))
 *
 * @param evaluate  Function pointer
 * @param n         Number of inputs
 * @return Newly allocated truth table
 */
TruthTable *nw_truth_table_create(bool (*evaluate)(const bool *, size_t), size_t n);

/**
 * @brief Evaluate a truth table at a specific input.
 * @param tt     Truth table
 * @param input  n-bit input
 * @return f(input)
 */
bool nw_truth_table_eval(const TruthTable *tt, const bool *input);

/**
 * @brief Free a truth table.
 * @param tt  Truth table to free
 */
void nw_truth_table_free(TruthTable *tt);

/**
 * @brief Check if a truth table represents a hard function.
 *
 * A function is hard if its truth table has no "simple structure"
 * detectable by examining sub-cubes or Fourier coefficients.
 *
 * This implements several heuristic tests:
 *   1. Balance check (Hamming weight near 2^{n-1})
 *   2. Autocorrelation check
 *   3. Sub-cube uniformity check
 *
 * @param tt      Truth table
 * @param thresh  Hardness threshold (circuit size)
 * @return true if the function appears hard
 */
bool nw_truth_table_is_hard(const TruthTable *tt, double thresh);

/**
 * @brief Compute the Fourier coefficient at subset S.
 *
 * \hat{f}(S) = (1/2^n) ∑_x f(x) · (-1)^{∑_{i∈S} x_i}
 *
 * L3: Mathematical Structure — Fourier analysis of boolean functions
 *
 * @param tt  Truth table
 * @param S   Subset indicator (size n)
 * @return Fourier coefficient ∈ [-1, 1]
 */
double nw_fourier_coefficient(const TruthTable *tt, const bool *S);

/**
 * @brief Compute the influence of variable i.
 *
 * Inf_i(f) = Pr_x[f(x) ≠ f(x^{⊕i})]
 * where x^{⊕i} is x with bit i flipped.
 *
 * L3: Mathematical Structure — Influence of variables
 *
 * @param tt  Truth table
 * @param i   Variable index
 * @return Influence ∈ [0, 1]
 */
double nw_variable_influence(const TruthTable *tt, size_t i);

/* ──────────────────────────────────────────────
 * L8: Natural Proofs Barrier
 * ────────────────────────────────────────────── */

/**
 * @brief Check if a given hardness result falls under the
 * natural proofs barrier.
 *
 * Razborov-Rudich (1997): No natural proof can show super-polynomial
 * circuit lower bounds unless no pseudorandom generators exist.
 *
 * A proof is "natural" if:
 *   1. Constructive: deciding the property is in P/poly
 *   2. Largeness: > 1/poly(n) fraction of functions have the property
 *   3. Usefulness: the property implies hardness
 *
 * L8: Advanced Topic — Natural proofs barrier
 *
 * @param hardness  Target hardness (circuit lower bound)
 * @param n         Input length
 * @return true if the natural proofs barrier applies
 */
bool nw_natural_proofs_barrier(double hardness, size_t n);

/**
 * @brief Compute the Razborov-Rudich threshold.
 *
 * The threshold circuit size above which natural proofs
 * cannot exist without breaking PRGs.
 *
 * @param n  Input length
 * @return Threshold circuit size
 */
double nw_razborov_rudich_threshold(size_t n);

/* ──────────────────────────────────────────────
 * L7: Derandomization via Hardness (IW97)
 * ────────────────────────────────────────────── */

/**
 * @brief Run Impagliazzo-Wigderson derandomization level analysis.
 *
 * IW97 showed three levels of derandomization from hardness:
 *   Level 0 (weak):    BPP ⊆ SUBEXP
 *   Level 1 (moderate): BPP ⊆ QP = DTIME(2^{polylog(n)})
 *   Level 2 (strong):  BPP = P
 *
 * L7: Application — BPP derandomization
 * L4: Fundamental Law — Impagliazzo-Wigderson Theorem
 *
 * @param level      Derandomization level (0, 1, or 2)
 * @param n          Problem size
 * @param seed_len   Output: required seed length for PRG
 * @return true if derandomization is achievable at this level
 */
bool nw_iw_derandomize(int level, size_t n, size_t *seed_len);

/**
 * @brief Compute IW97 sample size for BPP simulation.
 *
 * Given a BPP algorithm with error ε, how many samples from the
 * NW PRG are needed to simulate it deterministically?
 *
 * L7: Application — BPP simulation via NW PRG
 *
 * @param n        Input size
 * @param epsilon  Error parameter
 * @return Number of PRG samples needed
 */
size_t nw_bpp_simulation_samples(size_t n, double epsilon);

/* ──────────────────────────────────────────────
 * Utility: Hardness from truth table analysis
 * ────────────────────────────────────────────── */

/**
 * @brief Estimate the circuit complexity of a function from its
 * truth table using heuristic measures.
 *
 * Uses:
 *   - Sub-cube distribution entropy
 *   - Fourier concentration
 *   - Autocorrelation decay
 *
 * @param tt  Truth table
 * @return Estimated minimum circuit size
 */
double nw_estimate_circuit_complexity(const TruthTable *tt);

#endif /* NW_HARDNESS_H */
