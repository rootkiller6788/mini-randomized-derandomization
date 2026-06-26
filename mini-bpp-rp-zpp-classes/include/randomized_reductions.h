/******************************************************************************
 * randomized_reductions.h — Randomized Reductions & Interactive Proofs
 *
 * Covers randomized reductions between problems and the connection
 * between randomized complexity classes and interactive proofs.
 *
 * Key topics:
 *   - Randomized Karp/Levin reductions
 *   - RP-reductions, BPP-reductions
 *   - MA (Merlin-Arthur) and AM (Arthur-Merlin) classes
 *   - IP = PSPACE connection
 *   - Graph non-isomorphism in AM
 *
 * L1: Definitions — randomized reduction
 * L2: Core Concepts — completeness under randomized reductions
 * L4: Fundamental Laws — AM ⊆ Π₂^p, IP = PSPACE
 * L6: Canonical Problems — GNI (Graph Non-Isomorphism)
 * L8: Advanced Topics — PCP and hardness of approximation
 *
 * References:
 *   Arora-Barak Ch.8 (Interactive Proofs), Ch.18 (PCP)
 *   Goldwasser-Micali-Rackoff (1989)
 *   Babai (1985); Goldwasser-Sipser (1986)
 ******************************************************************************/

#ifndef RANDOMIZED_REDUCTIONS_H
#define RANDOMIZED_REDUCTIONS_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include "randomized_classes.h"

/* ================================================================
 * L1: Randomized Reduction Types
 * ================================================================ */

/** @brief A randomized Karp reduction from L₁ to L₂.
 *
 * A function f: Σ* → Σ* computable in randomized polynomial time
 * such that:
 *   x ∈ L₁ ⇒ Pr[f(x) ∈ L₂] ≥ 2/3
 *   x ∉ L₁ ⇒ Pr[f(x) ∉ L₂] ≥ 2/3  (BPP reduction)
 *
 * Or the one-sided variant (RP reduction):
 *   x ∈ L₁ ⇒ Pr[f(x) ∈ L₂] ≥ 1/2
 *   x ∉ L₁ ⇒ f(x) ∉ L₂ always
 */
typedef struct {
    size_t (*map)(const char *input, size_t input_len,
                  char *output, size_t output_cap,
                  RandomSource *rs);
    RandomizedClass reduction_type;  /**< BPP or RP reduction */
    size_t input_blowup;             /**< |f(x)| = poly(|x|) */
} RandomizedReduction;

/** @brief A language L with its decider and known complexity class. */
typedef struct {
    const char *name;
    bool (*decider)(const char *input, size_t len, RandomSource *rs);
    RandomizedClass class_type;
    size_t (*deterministic_decider)(const char*, size_t); /**< If in P */
} RandomizedLanguage;

/* ================================================================
 * L2: Randomized Reduction API
 * ================================================================ */

/**
 * @brief Compose two randomized reductions.
 *
 * If L₁ ≤_r L₂ via f₁ and L₂ ≤_r L₃ via f₂,
 * then L₁ ≤_r L₃ via f₂ ∘ f₁ (with amplified error).
 *
 * @param r1   First reduction (L₁ → L₂)
 * @param r2   Second reduction (L₂ → L₃)
 * @param k    Amplification parameter for composition
 * @return     Composed reduction (L₁ → L₃)
 */
RandomizedReduction compose_reductions(
    const RandomizedReduction *r1,
    const RandomizedReduction *r2,
    size_t k);

/**
 * @brief Amplify a randomized reduction to error ≤ 2^{-k}.
 *
 * Uses parallel repetition + majority vote.
 */
RandomizedReduction amplify_reduction(
    const RandomizedReduction *r, size_t k);

/**
 * @brief Verify a randomized reduction is correct on a test set.
 *
 * @param reduction  The reduction to verify
 * @param inputs     Test inputs
 * @param labels     True labels (1 = in language)
 * @param num_tests  Number of test cases
 * @param oracle     Membership oracle for target language
 * @param trials     Trials per input
 * @return           Fraction of inputs where reduction was correct
 */
double verify_reduction(
    const RandomizedReduction *reduction,
    const char **inputs, const bool *labels, size_t num_tests,
    bool (*oracle)(const char*, size_t),
    size_t trials);

/* ================================================================
 * L4: Arthur-Merlin (AM) Protocols
 * ================================================================ */

/**
 * @brief Arthur-Merlin protocol for Graph Non-Isomorphism (GNI).
 *
 * AM protocol:
 * 1. Arthur sends random permutation of G₀ or G₁ (coin flip)
 * 2. Merlin guesses which graph Arthur permuted
 * 3. If Merlin can distinguish, G₀ ≇ G₁
 *
 * If G₀ ≇ G₁: Merlin can always identify (soundness = 1)
 * If G₀ ≅ G₁: Merlin guesses correctly with prob 1/2
 *   → by parallel repetition, soundness ≤ 2^{-k}
 *
 * Reference: Goldwasser-Micali-Rackoff (1989); Arora-Barak §8.1
 *
 * @param adj0       Adjacency matrix of G₀ (n×n)
 * @param adj1       Adjacency matrix of G₁ (n×n)
 * @param n          Number of vertices
 * @param rounds     Number of parallel rounds
 * @param merlin_fn  Merlin's strategy (oracle/person)
 * @param rs         Randomness source (Arthur's coins)
 * @return           Arthur's decision: true = non-isomorphic
 */
bool am_graph_non_isomorphism(
    const int *adj0, const int *adj1, size_t n,
    size_t rounds,
    int (*merlin_fn)(const int *challenged_graph, size_t n),
    RandomSource *rs);

/**
 * @brief Generic AM[k] protocol simulator.
 *
 * An AM[k] protocol has k rounds: Arthur, Merlin, Arthur, Merlin, ...
 * (Arthur always moves first in AM, Merlin in MA).
 *
 * @param k             Number of rounds
 * @param input         Input string
 * @param input_len     Input length
 * @param arthur_fn     Arthur's message generation
 * @param merlin_fn     Merlin's response function
 * @param verifier_fn   Final verification
 * @param rs            Randomness source
 * @return              Verdict (true = accept)
 */
bool am_protocol_simulate(
    size_t k,
    const char *input, size_t input_len,
    char *(*arthur_fn)(const char *history, size_t hist_len, RandomSource *rs),
    char *(*merlin_fn)(const char *history, size_t hist_len),
    bool (*verifier_fn)(const char *history, size_t hist_len),
    RandomSource *rs);

/* ================================================================
 * L6: Canonical Problems with Randomized Algorithms
 * ================================================================ */

/**
 * @brief Randomized algorithm for Global Min-Cut (Karger's algorithm).
 *
 * Karger's contraction algorithm: repeatedly contract random edges
 * until 2 vertices remain. The probability of finding the min cut
 * is ≥ 1/C(n,2) = 2/(n(n-1)).
 *
 * By repeating O(n² log n) times, we find the min cut with high probability.
 * This is an RP-type algorithm (one-sided error: never outputs cut smaller
 * than true min cut).
 *
 * Reference: Karger (1993); Arora-Barak §7.2
 *
 * @param adj_matrix  n×n adjacency matrix (weights as integers)
 * @param n           Number of vertices
 * @param reps        Number of repetitions
 * @param rs          Randomness source
 * @return            Size of minimum cut found
 */
int karger_min_cut(const int *adj_matrix, size_t n, size_t reps,
                   RandomSource *rs);

/**
 * @brief Randomized 2-SAT algorithm (Papadimitriou's random walk).
 *
 * Given a 2-CNF formula with n variables and m clauses, find a
 * satisfying assignment or report unsatisfiable.
 *
 * Algorithm: start with random assignment, repeatedly pick an
 * unsatisfied clause and flip one of its variables randomly.
 * Expected time: O(n²) for satisfiable formulas.
 *
 * This is an RP algorithm (one-sided error on unsatisfiable instances).
 *
 * Reference: Papadimitriou (1991); Arora-Barak §7.2
 *
 * @param clauses     Array of 2m integers: (var_i, neg_i, var_j, neg_j)
 *                     where neg=1 means negated, var in [0,n-1]
 * @param n           Number of variables
 * @param m           Number of clauses
 * @param max_steps   Maximum random walk steps
 * @param rs          Randomness source
 * @param assignment  Output: satisfying assignment (if found)
 * @return            true if satisfying assignment found
 */
bool randomized_2sat(
    const int *clauses, size_t n, size_t m,
    size_t max_steps, RandomSource *rs,
    bool *assignment);

/**
 * @brief Randomized 3-SAT Schöning's algorithm.
 *
 * Repeat: pick random assignment, then do 3n random walk steps.
 * Success probability per trial: ≥ (3/4)^n · (some poly factor).
 * Overall: O((4/3)^n · poly(n)) time.
 *
 * Reference: Schöning (1999)
 *
 * @param formula  3-CNF formula as array of clauses
 * @param n        Number of variables
 * @param m        Number of clauses
 * @param rs       Randomness source
 * @param assign   Output: satisfying assignment (if found)
 * @return         true if satisfied
 */
bool schoning_3sat(
    const int *formula, size_t n, size_t m,
    RandomSource *rs, bool *assign);

/* ================================================================
 * L8: PCP Theorem Connection
 * ================================================================ */

/**
 * @brief PCP verifier simulation (simplified).
 *
 * The PCP Theorem: NP = PCP_{1/2,1}(O(log n), O(1))
 *
 * Every NP statement has a probabilistically checkable proof that
 * can be verified by reading only O(1) bits with error ≤ 1/2.
 *
 * This function simulates a simplified PCP verifier: query q bits
 * of proof, run a 3-query test, accept/reject.
 *
 * Reference: Arora-Barak Ch.18
 *
 * @param proof       The proof string π
 * @param proof_len   Length of proof
 * @param input       The input/claim being verified
 * @param input_len   Input length
 * @param queries     Output: indices of queried positions
 * @param num_queries Number of queries (typically 3-5)
 * @param rs          Randomness source
 * @return            Verdict
 */
bool pcp_verifier_simulate(
    const bool *proof, size_t proof_len,
    const char *input, size_t input_len,
    size_t *queries, size_t num_queries,
    RandomSource *rs);

/**
 * @brief Construct a PCP proof for 3SAT (simplified, theoretical).
 *
 * The full construction uses: polynomial encoding, low-degree tests,
 * sum-check protocol, and composition. This function implements a
 * simplified version for demonstration.
 *
 * @param formula   3-CNF formula
 * @param n         Number of variables
 * @param m         Number of clauses
 * @param proof_out Output proof buffer
 * @param proof_cap Capacity of proof buffer
 * @return          Actual proof length
 */
size_t pcp_proof_construct(
    const int *formula, size_t n, size_t m,
    bool *proof_out, size_t proof_cap);

#endif /* RANDOMIZED_REDUCTIONS_H */
