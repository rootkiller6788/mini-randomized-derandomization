/******************************************************************************
 * chernoff_bounds.h — Chernoff / Hoeffding Bounds & Concentration Inequalities
 *
 * Probability amplification in randomized algorithms relies on concentration
 * inequalities. These bounds quantify how quickly the average of independent
 * random variables converges to the expectation.
 *
 * L2: Core Concepts — tail bounds for error analysis
 * L3: Math Structures — probability distributions, moment generating functions
 * L4: Fundamental Laws — Chernoff-Hoeffding theorems
 *
 * References:
 *   Arora-Barak Appendix A.1, §7.1
 *   Mitzenmacher & Upfal, Probability and Computing (2005)
 *   Hoeffding (1963), Chernoff (1952)
 ******************************************************************************/

#ifndef CHERNOFF_BOUNDS_H
#define CHERNOFF_BOUNDS_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>

/* ================================================================
 * L2-L3: Core Types
 * ================================================================ */

/** @brief A Bernoulli trial: outcome 0 or 1 with probability p. */
typedef struct {
    double p;           /**< Probability of success (outcome 1) */
    uint64_t trials;    /**< Number of trials conducted */
    uint64_t successes; /**< Number of successes observed */
    double empirical_p; /**< Empirical success frequency */
} BernoulliExperiment;

/** @brief Generic bounded random variable for Hoeffding bounds. */
typedef struct {
    double lo;          /**< Lower bound of support [lo, hi] */
    double hi;          /**< Upper bound of support */
    double *samples;    /**< Array of observed samples */
    size_t n;           /**< Number of samples */
    double mean;        /**< Empirical mean */
    double variance;    /**< Empirical variance */
} BoundedSamples;

/* ================================================================
 * L4: Chernoff Bound — Multiplicative form
 * ================================================================ */

/**
 * @brief Chernoff upper tail bound for sum of independent Bernoulli trials.
 *
 * Let X = Σ X_i where X_i ~ Bernoulli(p_i) independent.
 * Let μ = E[X] = Σ p_i. Then for δ > 0:
 *
 *   Pr[X ≥ (1+δ)μ] ≤ exp(-μ·δ²/(2+δ))
 *
 * and for 0 < δ < 1:
 *
 *   Pr[X ≤ (1-δ)μ] ≤ exp(-μ·δ²/2)
 *
 * Reference: Mitzenmacher-Upfal Theorem 4.4, 4.5
 *
 * @param mu    Expected value μ = E[X]
 * @param delta Relative deviation δ
 * @return      Upper tail bound Pr[X ≥ (1+δ)μ]
 */
double chernoff_upper_tail(double mu, double delta);

/**
 * @brief Chernoff lower tail bound.
 *
 * @param mu    Expected value μ
 * @param delta Relative deviation δ ∈ (0,1)
 * @return      Lower tail bound Pr[X ≤ (1-δ)μ]
 */
double chernoff_lower_tail(double mu, double delta);

/**
 * @brief Two-sided Chernoff bound.
 *
 * Pr[|X - μ| ≥ δμ] ≤ 2 exp(-μ·δ²/3) for δ ∈ (0,1)
 *
 * @param mu    Expected value
 * @param delta Relative deviation
 * @return      Two-sided bound
 */
double chernoff_two_sided(double mu, double delta);

/**
 * @brief Compute the number of trials needed for BPP amplification.
 *
 * From error 1/3 to target error δ using majority vote, we need
 * m ≥ (3/δ²) · ln(2/δ) independent repetitions.
 *
 * More precisely, for Chernoff with μ = (2/3)m (on correct answer),
 * we need m such that Pr[majority wrong] ≤ exp(-m/18) ≤ δ
 * ⇒ m ≥ 18 ln(1/δ).
 *
 * Reference: Arora-Barak Lemma 7.9
 *
 * @param delta Target error probability
 * @return      Required number of independent trials
 */
size_t bpp_trials_for_error(double delta);

/**
 * @brief Compute trials needed for RP amplification.
 *
 * For RP (one-sided error), k trials give error ≤ (1/2)^k.
 * So k = ⌈log₂(1/δ)⌉ suffices.
 *
 * @param delta Target error probability
 * @return      Required number of independent trials
 */
size_t rp_trials_for_error(double delta);

/* ================================================================
 * L4: Hoeffding Bound — Additive form (no knowledge of p)
 * ================================================================ */

/**
 * @brief Hoeffding's inequality for bounded independent random variables.
 *
 * Let X₁,...,X_n be independent with a_i ≤ X_i ≤ b_i.
 * Let S_n = Σ X_i. Then for t > 0:
 *
 *   Pr[S_n - E[S_n] ≥ t] ≤ exp(-2t² / Σ(b_i - a_i)²)
 *
 * For the special case where all X_i ∈ [0,1]:
 *
 *   Pr[|(1/n)ΣX_i - E[X]| ≥ ε] ≤ 2 exp(-2nε²)
 *
 * Reference: Hoeffding (1963), Arora-Barak Lemma A.7
 *
 * @param n        Number of samples
 * @param epsilon  Deviation from mean
 * @return         Upper bound on Pr[|X̄ - μ| ≥ ε]
 */
double hoeffding_bound(size_t n, double epsilon);

/**
 * @brief Compute n needed for confidence interval of width 2ε
 * with confidence 1-δ using Hoeffding bound.
 *
 * n ≥ ln(2/δ) / (2ε²)
 *
 * @param epsilon Half-width of confidence interval
 * @param delta   Allowed failure probability
 * @return        Required sample size
 */
size_t hoeffding_sample_size(double epsilon, double delta);

/* ================================================================
 * L5: Empirical Concentration Analysis
 * ================================================================ */

/**
 * @brief Initialize a Bernoulli experiment tracker.
 */
BernoulliExperiment bernoulli_experiment_create(double p);

/**
 * @brief Update experiment with a new trial outcome.
 */
void bernoulli_experiment_update(BernoulliExperiment *exp, int outcome);

/**
 * @brief Check if observed deviation is statistically significant
 * at significance level α using Chernoff bound.
 *
 * @param exp    The experiment
 * @param alpha  Significance level
 * @return       true if deviation from expected p is significant
 */
bool bernoulli_experiment_significant(
    const BernoulliExperiment *exp, double alpha);

/**
 * @brief Compute exact binomial tail probability.
 *
 * Pr[Bin(n,p) ≥ k] = Σ_{i=k}^n C(n,i) p^i (1-p)^{n-i}
 *
 * Uses log-space computation to avoid overflow.
 *
 * @param n  Number of trials
 * @param k  Threshold
 * @param p  Success probability
 * @return   Exact tail probability
 */
double binomial_tail_exact(size_t n, size_t k, double p);

/**
 * @brief Compute exact binomial tail via log-space to handle large n.
 */
double binomial_tail_logspace(size_t n, size_t k, double p);

/**
 * @brief Compare Chernoff bound to exact binomial tail for given parameters.
 *
 * @param n     Number of trials
 * @param p     Success probability
 * @param delta Relative deviation
 * @param chernoff_bound Output: Chernoff upper bound
 * @param exact_value    Output: Exact binomial tail
 */
void chernoff_vs_exact(size_t n, double p, double delta,
                       double *chernoff_bound, double *exact_value);

/* ================================================================
 * L4: Advanced Concentration — Azuma-Hoeffding for Martingales
 * ================================================================ */

/**
 * @brief Azuma-Hoeffding inequality for martingale difference sequences.
 *
 * Let M_0,...,M_n be a martingale with |M_i - M_{i-1}| ≤ c_i.
 * Then for λ > 0:
 *
 *   Pr[M_n - M_0 ≥ λ] ≤ exp(-λ² / (2 Σ c_i²))
 *
 * Reference: Azuma (1967); Mitzenmacher-Upfal Theorem 12.4
 *
 * @param c_sq_sum  Σ c_i² — sum of squared differences bounds
 * @param lambda    Deviation threshold
 * @return          Upper bound on Pr[M_n - M_0 ≥ λ]
 */
double azuma_hoeffding_bound(double c_sq_sum, double lambda);

/**
 * @brief McDiarmid's inequality (bounded differences).
 *
 * If function f satisfies the bounded differences property with
 * constants c_i, then:
 *
 *   Pr[f(X) - E[f(X)] ≥ t] ≤ exp(-2t² / Σ c_i²)
 *
 * @param c_sq_sum  Σ c_i²
 * @param t         Deviation threshold
 * @return          Upper tail bound
 */
double mcdiarmid_bound(double c_sq_sum, double t);

/* ================================================================
 * L5: Application — Derandomization via Method of Conditional Probabilities
 * ================================================================ */

/**
 * @brief Method of conditional expectations for derandomization.
 *
 * Given a randomized algorithm whose success probability is proven
 * via Chernoff/Hoeffding, the method of conditional probabilities
 * deterministically finds a good random seed by iteratively fixing
 * bits while maintaining the conditional expectation above the threshold.
 *
 * This transforms a BPP algorithm into a P algorithm when the
 * conditional expectations are efficiently computable.
 *
 * Reference: Arora-Barak §7.3; Raghavan (1988)
 *
 * @param n             Number of random bits
 * @param objective_fn  Function to maximize (returns E[value | prefix])
 * @param threshold     Target value to exceed
 * @param seed_out      Output: deterministic "good" seed of length n
 * @return              true if a good seed was found
 */
bool conditional_expectation_derandomize(
    size_t n,
    double (*objective_fn)(const bool *prefix, size_t prefix_len),
    double threshold,
    bool *seed_out);

#endif /* CHERNOFF_BOUNDS_H */
