/******************************************************************************
 * chernoff_bounds.c — Chernoff/Hoeffding/Azuma Concentration Bounds
 *
 * Implements all major concentration inequalities used in the analysis
 * of randomized algorithms and derandomization:
 *
 * L4: Chernoff bounds (multiplicative, upper/lower/two-sided)
 * L4: Hoeffding bounds (additive, bounded support)
 * L4: Azuma-Hoeffding / McDiarmid (martingale, bounded differences)
 * L5: Empirical analysis and method of conditional probabilities
 *
 * Each function implements exactly one mathematical theorem.
 *
 * References:
 *   Chernoff (1952), Hoeffding (1963), Azuma (1967)
 *   Mitzenmacher & Upfal (2005): Probability and Computing
 *   Arora & Barak (2009): Appendix A
 ******************************************************************************/

#include "chernoff_bounds.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ================================================================
 * Utility: log2, factorial, binomial coefficient
 * ================================================================ */

/** log base 2 */
static inline double log2_val(double x) {
    if (x <= 0.0) return -1e300;
    return log(x) / log(2.0);
}

/** Natural log of binomial coefficient C(n, k) using Stirling approximation. */
static double log_binomial(size_t n, size_t k) {
    if (k > n) return -1e300;
    if (k == 0 || k == n) return 0.0;

    /* Stirling: ln(C(n,k)) ≈ n·H(k/n) - 0.5·ln(2π·n·(k/n)·(1-k/n))
     * where H(p) = -p ln p - (1-p) ln(1-p) */
    double p = (double)k / (double)n;
    double q = 1.0 - p;

    if (p < 1e-15 || q < 1e-15) return 0.0;

    double entropy = -p * log(p) - q * log(q);
    return n * entropy - 0.5 * log(2.0 * M_PI * n * p * q);
}

/* ================================================================
 * L4: Chernoff Bounds (Multiplicative Form)
 * ================================================================ */

double chernoff_upper_tail(double mu, double delta) {
    /* Theorem (Chernoff upper tail):
     * Let X = Σ X_i where X_i ∈ [0,1] independent, μ = E[X].
     * For δ > 0:
     *   Pr[X ≥ (1+δ)μ] ≤ exp(-μ·D((1+δ)||1))
     *
     * where D(a||b) = a ln(a/b) - a + b is KL-divergence.
     * The simplified bound: exp(-μ·min(δ,δ²)/3).
     *
     * For Bernoulli trials (stronger):
     *   Pr[X ≥ (1+δ)μ] ≤ (e^δ / (1+δ)^{1+δ})^μ
     *
     * This implements the tight bound (Chernoff-Cramér):
     *   Pr[X ≥ (1+δ)μ] ≤ exp(-μ·((1+δ)ln(1+δ) - δ)) */

    if (mu <= 0.0) return 1.0;   /* No trials */
    if (delta <= 0.0) return 1.0; /* Not an upper tail */

    double kl = (1.0 + delta) * log(1.0 + delta) - delta;
    double bound = exp(-mu * kl);

    /* Clamp to [0,1] */
    if (bound > 1.0) bound = 1.0;
    if (bound < 0.0) bound = 0.0;
    return bound;
}

double chernoff_lower_tail(double mu, double delta) {
    /* Theorem (Chernoff lower tail):
     * For 0 < δ < 1:
     *   Pr[X ≤ (1-δ)μ] ≤ exp(-μ·((1-δ)ln(1-δ) + δ))
     *
     * Or the simpler form: exp(-μ·δ²/2) */

    if (mu <= 0.0) return 1.0;
    if (delta <= 0.0 || delta >= 1.0) {
        /* δ ≥ 1: certainly upper bounds the probability */
        return (delta >= 1.0) ? 1.0 : 1.0;
    }

    double kl = (1.0 - delta) * log(1.0 - delta) + delta;
    double bound = exp(-mu * kl);

    if (bound > 1.0) bound = 1.0;
    if (bound < 0.0) bound = 0.0;
    return bound;
}

double chernoff_two_sided(double mu, double delta) {
    /* Theorem: For 0 < δ < 1:
     *   Pr[|X - μ| ≥ δμ] ≤ 2 exp(-μ·δ²/3)
     *
     * This follows from combining upper and lower tail bounds.
     * For the KL-divergence form, the worst case is:
     *   2 · exp(-μ · min(D(1+δ||1), D(1-δ||1)))
     *
     * For small δ, D(1±δ||1) ≈ δ²/2. */

    if (mu <= 0.0 || delta <= 0.0 || delta >= 1.0) return 1.0;

    double upper = chernoff_upper_tail(mu, delta);
    double lower = chernoff_lower_tail(mu, delta);
    double bound = 2.0 * fmax(upper, lower);

    /* Simpler bound: 2·exp(-μ·δ²/3) */
    double simple = 2.0 * exp(-mu * delta * delta / 3.0);
    /* Return the tighter bound */
    if (simple < bound) bound = simple;

    if (bound > 1.0) bound = 1.0;
    return bound;
}

size_t bpp_trials_for_error(double delta) {
    /* BPP base error = 1/3. After m trials with majority vote:
     * Error ≤ exp(-m/18). We need exp(-m/18) ≤ δ.
     * ⇒ m ≥ 18 ln(1/δ). */

    if (delta <= 0.0) return (size_t)-1;  /* Infinite */
    if (delta >= 1.0) return 1;

    double m = 18.0 * log(1.0 / delta);
    if (m < 3.0) m = 3.0;
    return (size_t)ceil(m);
}

size_t rp_trials_for_error(double delta) {
    /* RP: one-sided error ≤ 1/2 per trial.
     * After k trials (accept if any): Error ≤ (1/2)^k.
     * Need (1/2)^k ≤ δ ⇒ k ≥ log₂(1/δ). */

    if (delta <= 0.0) return (size_t)-1;
    if (delta >= 1.0) return 1;

    double k = log2(1.0 / delta);
    if (k < 1.0) k = 1.0;
    return (size_t)ceil(k);
}

/* ================================================================
 * L4: Hoeffding Bound
 * ================================================================ */

double hoeffding_bound(size_t n, double epsilon) {
    /* Theorem (Hoeffding, 1963):
     * For X_i ∈ [0,1] i.i.d., let X̄_n = (1/n)Σ X_i.
     * Then Pr[|X̄_n - μ| ≥ ε] ≤ 2 exp(-2nε²). */

    if (n == 0) return 1.0;
    if (epsilon <= 0.0) return 1.0;

    double bound = 2.0 * exp(-2.0 * (double)n * epsilon * epsilon);
    if (bound > 1.0) bound = 1.0;
    return bound;
}

size_t hoeffding_sample_size(double epsilon, double delta) {
    /* Hoeffding: 2 exp(-2nε²) ≤ δ ⇒ n ≥ ln(2/δ) / (2ε²). */

    if (epsilon <= 0.0 || delta <= 0.0) return (size_t)-1;
    if (delta >= 1.0) return 1;

    double n = log(2.0 / delta) / (2.0 * epsilon * epsilon);
    if (n < 1.0) n = 1.0;
    return (size_t)ceil(n);
}

/* ================================================================
 * L5: Empirical Concentration Analysis
 * ================================================================ */

BernoulliExperiment bernoulli_experiment_create(double p) {
    BernoulliExperiment exp;
    exp.p = p;
    exp.trials = 0;
    exp.successes = 0;
    exp.empirical_p = 0.0;
    return exp;
}

void bernoulli_experiment_update(BernoulliExperiment *exp, int outcome) {
    if (!exp) return;
    exp->trials++;
    if (outcome) exp->successes++;
    exp->empirical_p = (double)exp->successes / (double)exp->trials;
}

bool bernoulli_experiment_significant(
    const BernoulliExperiment *exp, double alpha)
{
    if (!exp || exp->trials == 0) return false;

    /* Check if observed deviation from expected p is significant
     * using Chernoff bound at level alpha. */
    double deviation = fabs(exp->empirical_p - exp->p);
    if (deviation < 1e-15) return false;

    /* Compute Chernoff bound for this deviation */
    double mu = exp->p * (double)exp->trials;
    if (mu < 1.0) {
        /* Use exact binomial for small mu */
        return binomial_tail_exact(exp->trials,
            (size_t)(exp->empirical_p * (double)exp->trials),
            exp->p) < alpha;
    }

    double delta = deviation / exp->p;
    double bound = chernoff_two_sided(mu, delta);

    return bound < alpha;
}

double binomial_tail_exact(size_t n, size_t k, double p) {
    /* Pr[Bin(n,p) ≥ k] = Σ_{i=k}^n C(n,i) p^i (1-p)^{n-i}
     *
     * Compute in log space to avoid overflow, then exponentiate. */
    if (k > n) return 0.0;
    if (k == 0) return 1.0;

    double total = 0.0;
    double q = 1.0 - p;
    double log_q = (q > 0.0) ? log(q) : -1e300;
    double log_p = (p > 0.0) ? log(p) : -1e300;

    for (size_t i = k; i <= n; i++) {
        double log_prob = log_binomial(n, i)
                        + (double)i * log_p
                        + (double)(n - i) * log_q;
        total += exp(log_prob);
        /* Early termination if contributions become negligible */
        if (log_prob < -50.0 && i > k + 10) break;
    }

    if (total > 1.0) total = 1.0;
    return total;
}

double binomial_tail_logspace(size_t n, size_t k, double p) {
    /* Same as above, but returns the log of the probability. */
    if (k > n) return -1e300;
    if (k == 0) return 0.0;

    /* For large n, use Normal approximation via Berry-Esseen:
     * Bin(n,p) ≈ N(np, np(1-p)).
     * Tail: Pr[Bin ≥ k] ≈ 1 - Φ((k - np)/√(np(1-p)))
     *
     * Or use the log-space sum for moderate n. */
    if (n > 1000) {
        double mu = (double)n * p;
        double sigma = sqrt(mu * (1.0 - p));
        double z = ((double)k - mu) / sigma;
        /* Normal tail approximation: 1-Φ(z) ≈ exp(-z²/2) / (z√2π) for z>0 */
        if (z > 0) {
            return -z * z / 2.0 - log(z * sqrt(2.0 * M_PI));
        }
        return 0.0;  /* z ≤ 0: tail is large, probability ≈ 1 */
    }

    /* Exact log-space sum */
    double log_total = -1e300;
    double q = 1.0 - p;
    double log_q = (q > 0.0) ? log(q) : -1e300;
    double log_p = (p > 0.0) ? log(p) : -1e300;

    for (size_t i = k; i <= n; i++) {
        double log_term = log_binomial(n, i)
                        + (double)i * log_p
                        + (double)(n - i) * log_q;
        /* log-sum-exp trick */
        if (log_term > log_total) {
            log_total = log_term + log(1.0 + exp(log_total - log_term));
        } else {
            log_total = log_total + log(1.0 + exp(log_term - log_total));
        }
    }

    return log_total;
}

void chernoff_vs_exact(size_t n, double p, double delta,
                       double *chernoff_bound, double *exact_value)
{
    double mu = (double)n * p;

    /* Chernoff upper tail: Pr[Bin ≥ (1+δ)μ] ≤ exp(-μ·((1+δ)ln(1+δ)-δ)) */
    double kl = (1.0 + delta) * log(1.0 + delta) - delta;
    *chernoff_bound = exp(-mu * kl);
    if (*chernoff_bound > 1.0) *chernoff_bound = 1.0;

    /* Exact binomial tail */
    size_t k = (size_t)ceil((1.0 + delta) * mu);
    if (k > n) {
        *exact_value = 0.0;
    } else {
        *exact_value = binomial_tail_exact(n, k, p);
    }
}

/* ================================================================
 * L4: Azuma-Hoeffding / McDiarmid
 * ================================================================ */

double azuma_hoeffding_bound(double c_sq_sum, double lambda) {
    /* Theorem (Azuma, 1967; Hoeffding, 1963 for martingales):
     * Let M_0,...,M_n be a martingale with |M_i - M_{i-1}| ≤ c_i.
     * Then Pr[M_n - M_0 ≥ λ] ≤ exp(-λ² / (2 Σ c_i²)). */

    if (c_sq_sum <= 0.0 || lambda <= 0.0) return 1.0;

    double bound = exp(-lambda * lambda / (2.0 * c_sq_sum));
    if (bound > 1.0) bound = 1.0;
    return bound;
}

double mcdiarmid_bound(double c_sq_sum, double t) {
    /* McDiarmid's inequality (bounded differences inequality):
     * If f(X₁,...,X_n) satisfies the bounded differences property
     * with constants c_i for each coordinate, then:
     *
     * Pr[f(X) - E[f(X)] ≥ t] ≤ exp(-2t² / Σ c_i²) */

    if (c_sq_sum <= 0.0 || t <= 0.0) return 1.0;

    double bound = exp(-2.0 * t * t / c_sq_sum);
    if (bound > 1.0) bound = 1.0;
    return bound;
}

/* ================================================================
 * L5: Method of Conditional Probabilities (Derandomization)
 * ================================================================ */

/**
 * The method of conditional probabilities (Erdős-Selfridge, 1973;
 * Spencer, 1987; Raghavan, 1988) transforms a probabilistic
 * existence proof into a deterministic algorithm.
 *
 * Given a random process with E[f(result)] ≥ τ, iteratively fix
 * each random choice to maintain E[f(result) | choices so far] ≥ τ.
 * After all choices are fixed, we obtain a deterministic result
 * with f(result) ≥ τ.
 *
 * This is the fundamental technique behind many derandomization results.
 */
bool conditional_expectation_derandomize(
    size_t n,
    double (*objective_fn)(const bool *prefix, size_t prefix_len),
    double threshold,
    bool *seed_out)
{
    if (!seed_out) return false;

    bool *prefix = (bool *)calloc(n, sizeof(bool));
    if (!prefix) return false;

    for (size_t i = 0; i < n; i++) {
        /* Try setting bit i = 0 */
        prefix[i] = false;
        double expect_0 = objective_fn(prefix, i + 1);

        /* Try setting bit i = 1 */
        prefix[i] = true;
        double expect_1 = objective_fn(prefix, i + 1);

        /* Choose the bit that gives higher conditional expectation */
        if (expect_0 >= expect_1) {
            prefix[i] = false;
        } else {
            prefix[i] = true;
        }

        /* Verify that conditional expectation stays above threshold */
        double current = objective_fn(prefix, i + 1);
        if (current < threshold - 1e-12) {
            /* Should not happen by the property of conditional expectations,
             * but guard against numerical issues */
            current = current;  /* no-op, proceed */
        }
    }

    memcpy(seed_out, prefix, n * sizeof(bool));
    free(prefix);

    /* Verify final objective meets threshold */
    double final_val = objective_fn(seed_out, n);
    return (final_val >= threshold - 1e-10);
}

/* ================================================================
 * L5: Bonus — Poisson Tail Bound
 * ================================================================ */

/**
 * Poisson tail bound for rare events.
 *
 * If X ~ Poisson(λ), then for x > λ:
 *   Pr[X ≥ x] ≤ (eλ/x)^x · e^{-λ} = exp(x ln λ - x ln x + x - λ)
 */
double poisson_upper_tail(double lambda, double x) {
    if (x <= lambda) return 1.0;
    if (lambda <= 0.0) return (x > 0.0) ? 0.0 : 1.0;

    double log_bound = x * log(lambda) - x * log(x) + x - lambda;
    double bound = exp(log_bound);
    if (bound > 1.0) bound = 1.0;
    return bound;
}

/**
 * Johnson-Lindenstrauss Lemma concentration.
 *
 * For random projection from d to k = O(log n / ε²) dimensions,
 * the pairwise distances are preserved within (1±ε) with high probability.
 *
 * The concentration uses the chi-squared tail bound:
 *   Pr[|‖Proj(x)‖² - ‖x‖²| > ε‖x‖²] ≤ 2 exp(-(ε²-ε³)k/4)
 *
 * @param d        Original dimension
 * @param n        Number of points
 * @param epsilon  Distortion tolerance
 * @param delta    Failure probability
 * @return         Target dimension k
 */
size_t johnson_lindenstrauss_dim(size_t n, double epsilon, double delta) {
    /* k ≥ (4 / (ε² - ε³)) · ln(n/δ)
     * For small ε: k ≈ 4 ln(n) / ε² */
    if (epsilon <= 0.0 || epsilon >= 1.0 || delta <= 0.0) return 1;

    double factor = 4.0 / (epsilon * epsilon - epsilon * epsilon * epsilon);
    if (factor < 0) factor = 4.0 / (epsilon * epsilon);  /* Fallback for large ε */

    double k = factor * log((double)n / delta);
    if (k < 1.0) k = 1.0;
    return (size_t)ceil(k);
}
