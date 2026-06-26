/******************************************************************************
 * derandomization_methods.h — Derandomization & the BPP = P Question
 *
 * Covers key approaches to derandomization:
 *   1. Method of conditional probabilities (deterministic search)
 *   2. Enumeration of all random seeds (brute force, exponential)
 *   3. Pseudorandom generators (if E requires exponential circuits)
 *   4. Pairwise independence and limited independence
 *
 * L4: Fundamental Laws — BPP ⊆ P/poly (Adleman), BPP ⊆ Σ₂^p (Sipser-Gács)
 * L5: Algorithms/Methods — derandomization techniques
 * L7: Applications — cryptography, algorithm design
 * L8: Advanced Topics — hardness vs randomness paradigm
 *
 * References:
 *   Arora-Barak Ch.7, Ch.20 (Derandomization)
 *   Impagliazzo & Wigderson (1997)
 *   Vadhan, Pseudorandomness (2012)
 ******************************************************************************/

#ifndef DERANDOMIZATION_METHODS_H
#define DERANDOMIZATION_METHODS_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include "randomized_classes.h"
#include "chernoff_bounds.h"

/* ================================================================
 * L4: Method of Conditional Probabilities
 * ================================================================ */

/**
 * @brief Pessimistic estimator for derandomizing MAX-CUT
 * via the method of conditional probabilities.
 *
 * The randomized MAX-CUT algorithm assigns each vertex to left/right
 * uniformly at random, achieving E[|cut|] ≥ |E|/2.
 *
 * The pessimistic estimator maintains a lower bound on the expected
 * cut size given the partial assignment, allowing deterministic
 * selection of each vertex's side while preserving the guarantee.
 *
 * Reference: Raghavan (1988); Arora-Barak §7.3
 *
 * @param adj_matrix  n×n adjacency matrix (undirected, symmetric)
 * @param n           Number of vertices
 * @param assignments Output: binary assignment for each vertex
 * @return            Guaranteed cut size
 */
size_t max_cut_derandomized(
    const int *adj_matrix, size_t n,
    bool *assignments);

/**
 * @brief Derandomize set balancing via conditional probabilities.
 *
 * Given m subsets of {1..n}, assign ±1 to each element to minimize
 * the maximum discrepancy. The randomized assignment (random ±1)
 * gives discrepancy O(√(n log m)).
 *
 * Conditional probabilities gives a deterministic algorithm achieving
 * the same bound when m = poly(n).
 *
 * @param sets        m×n incidence matrix (1 if element j ∈ set i)
 * @param m           Number of sets
 * @param n           Number of elements
 * @param assignment  Output: ±1 for each element
 * @return            Maximum discrepancy achieved
 */
size_t set_balancing_derandomized(
    const int *sets, size_t m, size_t n,
    int *assignment);

/* ================================================================
 * L4: Pairwise Independent Hash Functions
 * ================================================================ */

/**
 * @brief Pairwise independent hash family H: [N] → [M].
 *
 * A family H = {h_a,b : x ↦ (ax + b) mod p mod M} where p ≥ N is prime,
 * a ∈ [1, p-1], b ∈ [0, p-1].
 *
 * For any x₁ ≠ x₂, Pr_h[h(x₁)=y₁ ∧ h(x₂)=y₂] = 1/M².
 *
 * This reduces randomness from O(N log M) to O(log N + log M) bits.
 *
 * Reference: Carter-Wegman (1979); Arora-Barak §7.3.2
 */
typedef struct {
    uint64_t a;         /**< Multiplier */
    uint64_t b;         /**< Additive shift */
    uint64_t p;         /**< Prime modulus ≥ N */
    uint64_t M;         /**< Range size */
} PairwiseHash;

/**
 * @brief Create a random hash function from the pairwise independent family.
 *
 * @param N    Universe size (will find prime ≥ N)
 * @param M    Range size
 * @param rs   Randomness source (only uses 2 log N random bits!)
 * @return     A pairwise independent hash function
 */
PairwiseHash pairwise_hash_create(uint64_t N, uint64_t M, RandomSource *rs);

/**
 * @brief Evaluate the hash: h(x) = (a*x + b) mod p mod M
 */
uint64_t pairwise_hash_eval(const PairwiseHash *h, uint64_t x);

/**
 * @brief Verify pairwise independence: for given x₁≠x₂, y₁, y₂,
 * count how many keys (a,b) map accordingly.
 */
size_t pairwise_hash_verify(
    uint64_t p, uint64_t M,
    uint64_t x1, uint64_t x2,
    uint64_t y1, uint64_t y2);

/* ================================================================
 * L4: k-wise Independent Hash Functions
 * ================================================================ */

/**
 * @brief k-wise independent hash family via random polynomials over GF(p).
 *
 * h(x) = (Σ_{i=0}^{k-1} a_i x^i) mod p mod M
 *
 * Uses only k log p random bits instead of N log M.
 *
 * Reference: Joffe (1974); Wegman-Carter (1981)
 */
typedef struct {
    uint64_t *coeffs;   /**< Polynomial coefficients a_0,...,a_{k-1} */
    size_t k;           /**< Degree + 1 = k-wise independence */
    uint64_t p;         /**< Prime modulus */
    uint64_t M;         /**< Range size */
} KWiseHash;

/**
 * @brief Create a k-wise independent hash function.
 *
 * @param N    Universe size
 * @param M    Range size
 * @param k    Independence parameter
 * @param rs   Randomness source
 * @return     k-wise independent hash function
 */
KWiseHash kwise_hash_create(uint64_t N, uint64_t M, size_t k, RandomSource *rs);

/**
 * @brief Evaluate k-wise independent hash.
 */
uint64_t kwise_hash_eval(const KWiseHash *h, uint64_t x);

/**
 * @brief Free resources associated with the hash.
 */
void kwise_hash_destroy(KWiseHash *h);

/* ================================================================
 * L5: Deterministic Simulation of BPP (Exponential)
 * ================================================================ */

/**
 * @brief Brute-force derandomization: iterate over all 2^m random strings.
 *
 * This yields a deterministic algorithm with running time 2^m · T(n),
 * showing BPP ⊆ EXP trivially.
 *
 * For an RP machine, we need only search for an accepting random string;
 * for a BPP machine, we take majority over all random strings.
 *
 * @param input       Input string
 * @param input_len   Input length
 * @param decider     The probabilistic decision function
 * @param m           Number of random bits used by the decider
 * @param time_bound  Running time bound per trial
 * @param class_type  BPP (majority vote) or RP (exists accept)
 * @return            Deterministic decision
 */
bool brute_force_derandomize(
    const char *input, size_t input_len,
    bool (*decider)(const char*, size_t, RandomSource*),
    size_t m, size_t time_bound, RandomizedClass class_type);

/* ================================================================
 * L7: Application — Randomized Quicksort Analysis
 * ================================================================ */

/**
 * @brief Analyze randomized quicksort as a BPP/ZPP-type algorithm.
 *
 * Randomized quicksort with random pivot selection has:
 * - Expected running time O(n log n) → like ZPP (expected poly time)
 * - Tail bound: Pr[T > c·n log n] ≤ n^{-Θ(c)} → concentration
 *
 * This function models the probabilistic analysis.
 *
 * @param n          Number of elements
 * @param num_trials Number of simulation runs
 * @param mean       Output: empirical mean comparisons
 * @param stddev     Output: empirical standard deviation
 */
void randomized_quicksort_analysis(
    size_t n, size_t num_trials,
    double *mean, double *stddev);

/**
 * @brief Compare random vs median-of-3 pivot selection in quicksort.
 *
 * @return Ratio of median-of-3 comparisons to random pivot comparisons
 */
double quicksort_pivot_comparison(size_t n, size_t trials);

/* ================================================================
 * L8: Advanced — Hardness vs Randomness
 * ================================================================ */

/**
 * @brief The NW (Nisan-Wigderson) generator concept:
 * If there exists a function f in E = DTIME(2^{O(n)}) with circuit
 * complexity 2^{Ω(n)}, then BPP = P.
 *
 * This function encodes the core insight: hardness ⇒ derandomization.
 *
 * @param circuit_size_lower_bound  The assumed circuit lower bound
 * @param n                         Input length parameter
 * @return                          Required PRG seed length
 */
size_t nw_generator_seed_length(
    double circuit_size_lower_bound, size_t n);

/**
 * @brief Impagliazzo's five worlds: characterize the derandomization
 * implications of each world for P vs BPP.
 *
 * 1. Algorithmica:  P = BPP (or P = NP)
 * 2. Heuristica:    NP hard on average, BPP = P
 * 3. Pessiland:     OWF exist, but no public-key crypto
 * 4. Minicrypt:     OWF exist, BPP ⊂ EXP
 * 5. Cryptomania:   Public-key crypto exists
 *
 * This function maps a "world" to the truth of BPP = P in that world.
 *
 * @param world  World index (0-4)
 * @return       Whether BPP = P in this world
 */
typedef enum {
    IMP_ALGORITHMICA = 0,
    IMP_HEURISTICA   = 1,
    IMP_PESSILAND    = 2,
    IMP_MINICRYPT    = 3,
    IMP_CRYPTOMANIA  = 4,
} ImpagliazzoWorld;

const char *impagliazzo_world_name(ImpagliazzoWorld w);
bool bpp_equals_p_in_world(ImpagliazzoWorld w);
const char *impagliazzo_world_description(ImpagliazzoWorld w);

#endif /* DERANDOMIZATION_METHODS_H */
