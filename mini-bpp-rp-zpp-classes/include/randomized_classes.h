/******************************************************************************
 * randomized_classes.h — Core Randomized Complexity Classes
 *
 * Defines the fundamental probabilistic complexity classes:
 *   BPP (Bounded-error Probabilistic Polynomial time)
 *   RP  (Randomized Polynomial time — one-sided error)
 *   ZPP (Zero-error Probabilistic Polynomial time)
 *   co-RP, co-BPP
 *
 * Theorem sources: Arora & Barak (2009), Ch.7; Sipser (2013), Ch.10
 *
 * L1: Definitions — struct definitions for randomized Turing machines
 * L2: Core Concepts — error probability bounds, amplification
 * L3: Math Structures — probabilistic TM formalization
 ******************************************************************************/

#ifndef RANDOMIZED_CLASSES_H
#define RANDOMIZED_CLASSES_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>

/* ================================================================
 * L1: Core Data Structures
 * ================================================================ */

/** @brief Randomness source abstraction for probabilistic Turing machines.
 *
 * Models a fair coin-tossing mechanism with independent tosses.
 * Each call to randbit() returns 0 or 1 with probability 1/2 each,
 * independent of all previous calls.
 *
 * Reference: Arora-Barak §7.1, Definition 7.1
 */
typedef struct {
    uint64_t    seed;
    uint64_t    state;
    uint64_t    tosses_used;
    bool        exhausted;
} RandomSource;

/** @brief Probabilistic Turing Machine with explicit randomness tape. */
typedef struct {
    size_t      input_length;
    size_t      random_tape_len;
    double      error_prob;
    bool        (*decider)(const char *input, size_t len, RandomSource *r);
    size_t      time_bound;
    size_t      steps_taken;
} ProbabilisticTM;

/** @brief Outcome of a single probabilistic computation run. */
typedef enum {
    PTM_CORRECT_ACCEPT = 0,
    PTM_CORRECT_REJECT = 1,
    PTM_FALSE_REJECT   = 2,
    PTM_FALSE_ACCEPT   = 3,
} PTMOutcome;

/** @brief Classification of a language into randomized complexity classes. */
typedef enum {
    CLASS_BPP = 0,
    CLASS_RP  = 1,
    CLASS_CO_RP = 2,
    CLASS_ZPP = 3,
    CLASS_P   = 4,
    CLASS_UNKNOWN = 5,
} RandomizedClass;

/** @brief Result of a language membership test with confidence bounds. */
typedef struct {
    bool        accepted;
    double      confidence;
    size_t      trials;
    size_t      accept_count;
    size_t      reject_count;
    size_t      unknown_count;
} RandomizedDecision;

/* ================================================================
 * L2: Core API — Probability Amplification & Error Reduction
 * ================================================================ */

RandomSource random_source_init(uint64_t seed);
int random_source_bit(RandomSource *rs);
uint64_t random_source_uniform(RandomSource *rs, uint64_t n);
void random_source_reset(RandomSource *rs, uint64_t new_seed);

RandomizedDecision bpp_amplify(
    const char *input, size_t input_len,
    bool (*base_decider)(const char*, size_t, RandomSource*),
    size_t k, size_t time_bound);

RandomizedDecision rp_amplify(
    const char *input, size_t input_len,
    bool (*rp_decider)(const char*, size_t, RandomSource*),
    size_t k, size_t time_bound);

RandomizedDecision zpp_simulate(
    const char *input, size_t input_len,
    int (*zpp_decider)(const char*, size_t, RandomSource*),
    size_t max_attempts, size_t time_bound);

/* ================================================================
 * L3: Mathematical Structures — Formal PTM Simulation
 * ================================================================ */

ProbabilisticTM ptm_create(
    bool (*decider)(const char*, size_t, RandomSource*),
    size_t random_len, size_t time_bound, double error_prob);

PTMOutcome ptm_run_once(
    ProbabilisticTM *ptm, const char *input, size_t input_len,
    RandomSource *rs, bool oracle(const char*, size_t));

double ptm_estimate_accept_prob(
    ProbabilisticTM *ptm, const char *input, size_t input_len,
    size_t num_samples, double *conf_radius);

/* ================================================================
 * L4: Fundamental Laws — Class Relationships
 * ================================================================ */

bool verify_class_hierarchy(void);
void classify_machine(double false_pos_rate, double false_neg_rate,
                      double confidence, bool class_results[4]);
size_t adleman_theorem_advice_size(size_t n, double err_bound);
bool bpp_in_sigma2(void);

/* ================================================================
 * L5: Algorithms — Concrete BPP/RP/ZPP Implementations
 * ================================================================ */

bool pit_schwartz_zippel(const int *coefficients, size_t num_vars,
    size_t max_degree, size_t field_size, RandomSource *rs);

bool rp_perfect_matching(const int *adj_matrix, size_t n, RandomSource *rs);

int zpp_primality_test(uint64_t n, RandomSource *rs, bool deterministic);

bool freivalds_verify(const double *A, const double *B,
    const double *C, size_t n, RandomSource *rs);

bool freivalds_verify_amplified(const double *A, const double *B,
    const double *C, size_t n, size_t k, RandomSource *rs);

double approximate_count_sat(const char *formula, size_t len,
    double epsilon, double delta, RandomSource *rs);

#endif /* RANDOMIZED_CLASSES_H */
