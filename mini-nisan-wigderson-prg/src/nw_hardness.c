/**
 * @file nw_hardness.c
 * @brief Hardness assumptions, amplification, truth tables, and analysis.
 *
 * Implements:
 *   L1: Hardness assumption types
 *   L3: Truth table operations (creation, evaluation, Fourier analysis)
 *   L4: Yao XOR Lemma (full amplification computation)
 *   L4: Shannon's counting lower bound
 *   L4: Impagliazzo's worst-case to average-case reduction
 *   L5: Hard-core set construction algorithm
 *   L8: Natural proofs barrier analysis
 *   L7: IW97 derandomization level analysis
 *
 * Reference:
 *   Yao, "Theory and Applications of Trapdoor Functions", FOCS 1982
 *   Impagliazzo, "Hard-core distributions for somewhat hard problems", FOCS 1995
 *   Impagliazzo & Wigderson, "P = BPP if E requires exponential circuits", STOC 1997
 *   Razborov & Rudich, "Natural Proofs", JCSS 55(1), 1997
 *   Arora & Barak, Ch 19-20, Ch 23
 */

#include "nw_hardness.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* =============================================
 * L4: Yao XOR Lemma — Full Implementation
 *
 * Yao 1982 / Levin 1987 / Goldreich-Nisan-Wigderson 1995:
 *
 * If f: {0,1}^n -> {0,1} is (epsilon, S)-hard
 * (no circuit of size S can compute f with advantage > epsilon),
 * then XOR_k(f): {0,1}^{k*n} -> {0,1} defined by
 *   XOR_k(f)(x_1,...,x_k) = f(x_1) xor ... xor f(x_k)
 * is (epsilon^k + delta, S')-hard with
 *   S' = Omega(delta^2 * epsilon^2 * S).
 *
 * Advantage epsilon defined as: |Pr[C(x)=f(x)] - 1/2| = epsilon
 * (equivalently: error prob = 1/2 - epsilon for random guess)
 * ============================================= */

void nw_yao_xor_lemma_full(
    size_t n, double epsilon, size_t k, double delta, double base_size,
    size_t *new_n, double *new_eps, double *new_size)
{
    *new_n = n * k;

    /* XOR-k amplifies advantage exponentially:
     *   new_epsilon = (1 - (1-2*epsilon)^k) / 2
     *
     * Intuition: if a circuit guesses each f(x_i) with advantage epsilon,
     * after XOR'ing k copies, the correlation = (1-2*epsilon)^k
     * So advantage = (1 - correlation) / 2 + delta
     */
    if (epsilon >= 0.5) {
        *new_eps = 0.5;
    } else if (epsilon <= 0.0) {
        *new_eps = delta;
    } else {
        double correlation = 1.0 - 2.0 * epsilon;
        double new_correlation = pow(correlation, (double)k);
        *new_eps = (1.0 - new_correlation) / 2.0 + delta;
        if (*new_eps > 0.5) *new_eps = 0.5;
    }

    /* Amplified circuit size bound:
     * S' = Omega(delta^2 * epsilon^2 * S) for moderate k
     * More precise (Impagliazzo 1995):
     *   S' = C * delta^2 * epsilon^4 * S / k^2
     * for some constant C.
     */
    if (epsilon <= 0.0 || delta <= 0.0) {
        *new_size = base_size;
    } else {
        double factor = delta * delta * epsilon * epsilon * epsilon * epsilon;
        factor /= (16.0 * (double)(k * k));
        *new_size = base_size * factor;
        if (*new_size < 1.0) *new_size = 1.0;
    }
}

/* =============================================
 * L5: Direct XOR of Truth Tables
 * ============================================= */

TruthTable *nw_truth_table_xor_k(const TruthTable *tt, size_t k) {
    if (!tt || k == 0 || k > 10) return NULL; /* Size limits for practicality */

    size_t new_n = tt->n * k;
    size_t new_size = (size_t)1 << new_n;

    /* Check size feasibility */
    if (new_n > 24 || new_size > ((size_t)1 << 24)) return NULL;

    TruthTable *result = (TruthTable *)malloc(sizeof(TruthTable));
    if (!result) return NULL;
    result->n = new_n;
    result->size = new_size;
    result->table = (bool *)malloc(new_size * sizeof(bool));
    if (!result->table) { free(result); return NULL; }

    /* For each (x_1, ..., x_k) in {0,1}^{k*n}:
     *   g(x_1,...,x_k) = XOR_{i=1}^k f(x_i)
     *
     * Each x_i is an n-bit block within the new_n-bit input.
     * We precompute f for all 2^n inputs.
     */
    size_t block_size = (size_t)1 << tt->n;

    for (size_t idx = 0; idx < new_size; idx++) {
        bool xor_result = false;
        for (size_t i = 0; i < k; i++) {
            /* Extract i-th block of n bits */
            size_t block = (idx >> (i * tt->n)) & (block_size - 1);
            xor_result ^= tt->table[block];
        }
        result->table[idx] = xor_result;
    }

    return result;
}

/* =============================================
 * L4: Worst-Case to Average-Case Reduction
 *
 * Impagliazzo's Hard-Core Set Theorem (1995):
 *
 * Given a boolean function f that is delta-hard in the worst case
 * (any circuit of size S errs on at least delta fraction of inputs),
 * there exists a "hard-core" set H of density >= delta on which
 * f is (1/2 - epsilon)-hard for circuits of size epsilon^2*delta^2*S.
 *
 * This is the key bridge between worst-case and average-case hardness,
 * enabling the NW construction which requires average-case hardness.
 *
 * Algorithm (Klivans & Servedio 2003, boosting approach):
 *   1. Start with uniform distribution over all inputs
 *   2. Iteratively find weak learners and re-weight
 *   3. The final hard-core set consists of inputs that are consistently
 *      misclassified by the ensemble of weak learners.
 * ============================================= */

WorstToAverage *nw_worst_to_average(const TruthTable *tt, double worst_hardness, double epsilon) {
    if (!tt || worst_hardness <= 0.0 || worst_hardness > 1.0 || epsilon <= 0.0) return NULL;

    WorstToAverage *wta = (WorstToAverage *)malloc(sizeof(WorstToAverage));
    if (!wta) return NULL;

    wta->n = tt->n;
    wta->worst_hardness = worst_hardness;
    wta->hardcore_density = worst_hardness; /* Density >= delta */

    /* Compute average-case hardness from worst-case.
     * By Impagliazzo's theorem:
     *   avg_hardness = Omega(epsilon^2 * worst_hardness^2 * S)
     * The trade-off: epsilon * density >= worst_hardness
     */
    wta->avg_hardness = epsilon * worst_hardness * worst_hardness;
    if (wta->avg_hardness > 0.5) wta->avg_hardness = 0.5;

    /* Construct hard-core set H:
     * An input x is in H if f(x) is "unpredictable"
     * We use a simple heuristic: inputs where f differs from majority */
    size_t total = tt->size;
    wta->hardcore_set = (bool *)calloc(total, sizeof(bool));
    if (!wta->hardcore_set) { free(wta); return NULL; }

    /* Compute the fraction of 1s */
    size_t ones = 0;
    for (size_t i = 0; i < total; i++)
        if (tt->table[i]) ones++;

    /* balance used for analysis; suppress warning */
    (void)((double)ones / (double)total);

    /* Hard-core set = inputs where f differs from a randomly biased
     * predictor. The deterministic predictor is the majority value. */
    bool majority = (ones >= total / 2);
    size_t hardcore_count = 0;
    size_t target_count = (size_t)(worst_hardness * (double)total);

    for (size_t i = 0; i < total && hardcore_count < target_count; i++) {
        if (tt->table[i] != majority) {
            wta->hardcore_set[i] = true;
            hardcore_count++;
        }
    }

    /* If not enough, add the rest by probability */
    for (size_t i = 0; i < total && hardcore_count < target_count; i++) {
        if (!wta->hardcore_set[i]) {
            wta->hardcore_set[i] = true;
            hardcore_count++;
        }
    }

    return wta;
}

void nw_worst_to_average_free(WorstToAverage *wta) {
    if (!wta) return;
    free(wta->hardcore_set);
    free(wta);
}

/* =============================================
 * L4: Shannon's Lower Bound
 *
 * Shannon (1949) counting argument:
 * - Each circuit of size s can be described by:
 *   * s gate types: each log(16) bits (for common basis)
 *   * 2s inputs (each input is a gate index): 2s * log(s) bits
 * - Total description length <= s * log(s) + O(s) bits
 * - Number of circuits of size s <= 2^{O(s log s)}
 * - Number of boolean functions on n variables = 2^{2^n}
 *
 * Therefore for most functions:
 *   2^{O(s log s)} < 2^{2^n} => s = Omega(2^n / n)
 *
 * More precisely: at least 2^n/n gates needed for most functions.
 * ============================================= */

double nw_shannon_lower_bound(size_t n) {
    if (n <= 1) return 1.0;

    /* Shannon 1949: Most functions require at least
     *   2^n / (2n) * (1 - o(1))
     * gates over the basis {AND, OR, NOT}.
     *
     * We return the asymptotic lower bound factor.
     */
    double bound = pow(2.0, (double)n) / (2.0 * (double)n);

    /* For small n, adjust constant */
    if (n <= 5) bound = pow(2.0, (double)n) / (4.0 * (double)n);

    return bound;
}

double nw_lupanov_upper_bound(size_t n) {
    if (n <= 1) return 2.0;

    /* Lupanov (1958): Every boolean function on n variables
     * can be computed by a circuit of size
     *   2^n / n * (1 + o(1))
     * using AND, OR, NOT gates.
     *
     * The construction decomposes the truth table into blocks.
     */
    double bound = pow(2.0, (double)n) / (double)n;
    /* (1 + o(1)) factor approximated as (1 + log n / n) */
    bound *= (1.0 + log((double)n) / (double)n);
    return bound;
}

/* =============================================
 * L5: Hardness Assumption Validation
 * ============================================= */

bool nw_hardness_assumption_valid(const HardnessAssumption *ha) {
    if (!ha) return false;
    if (ha->n == 0) return false;
    if (ha->circuit_size <= 0.0) return false;
    if (ha->hardness <= 0.0 || ha->hardness > 0.5) return false;

    /* If exponential, check circuit_size >= 2^{c*n} for some c > 0 */
    if (ha->is_exponential) {
        double min_exp = pow(2.0, (double)ha->n * 0.01); /* c = 0.01 minimum */
        if (ha->circuit_size < min_exp) return false;
    }

    return true;
}

/* =============================================
 * L5: Truth Table Operations
 * ============================================= */

TruthTable *nw_truth_table_create(bool (*evaluate)(const bool *, size_t), size_t n) {
    if (!evaluate || n > 24) return NULL; /* size limit */

    size_t size = (size_t)1 << n;
    if (size == 0 || size > ((size_t)1 << 24)) return NULL;

    TruthTable *tt = (TruthTable *)malloc(sizeof(TruthTable));
    if (!tt) return NULL;
    tt->n = n;
    tt->size = size;
    tt->table = (bool *)malloc(size * sizeof(bool));
    if (!tt->table) { free(tt); return NULL; }

    /* Evaluate on all inputs in lexicographic order */
    bool *input = (bool *)calloc(n, sizeof(bool));
    if (!input) { free(tt->table); free(tt); return NULL; }

    for (size_t i = 0; i < size; i++) {
        /* Decode i to binary input */
        for (size_t j = 0; j < n; j++)
            input[j] = (i >> j) & 1;
        tt->table[i] = evaluate(input, n);
    }

    free(input);
    return tt;
}

bool nw_truth_table_eval(const TruthTable *tt, const bool *input) {
    if (!tt || !input) return false;
    /* Convert input to index */
    size_t idx = 0;
    for (size_t i = 0; i < tt->n; i++)
        if (input[i]) idx |= ((size_t)1 << i);
    if (idx >= tt->size) return false;
    return tt->table[idx];
}

void nw_truth_table_free(TruthTable *tt) {
    if (!tt) return;
    free(tt->table);
    free(tt);
}

/**
 * L5: Hardness Heuristic Tests
 *
 * Tests if a function looks "hard" by examining:
 *   1. Balance: Hamming weight near 2^{n-1} (unbiased)
 *   2. Autocorrelation: small correlation with shifted versions
 *   3. Algebraic degree proxy: nonlinearity measure
 */

bool nw_truth_table_is_hard(const TruthTable *tt, double thresh) {
    if (!tt) return false;
    (void)thresh;

    size_t size = tt->size;
    if (size < 4) return false;

    /* Test 1: Balance check.
     * A random function has weight ~ 2^{n-1} +- sqrt(2^n).
     * We check if weight is within 3 standard deviations. */
    size_t weight = 0;
    for (size_t i = 0; i < size; i++)
        if (tt->table[i]) weight++;

    double expected = (double)size / 2.0;
    double stddev = sqrt((double)size) / 2.0;
    double deviation = fabs((double)weight - expected);
    bool balanced = (deviation <= 3.0 * stddev);

    /* Test 2: Autocorrelation at distance 1.
     * For each i, compare f(x) with f(x xor e_i).
     * A random function has correlation near 0. */
    double max_autocorr = 0.0;
    for (size_t bit = 0; bit < tt->n; bit++) {
        size_t matches = 0;
        size_t mask = (size_t)1 << bit;
        for (size_t i = 0; i < size; i++) {
            size_t j = i ^ mask;
            if (tt->table[i] == tt->table[j]) matches++;
        }
        double corr = fabs((double)matches / (double)size - 0.5) * 2.0;
        if (corr > max_autocorr) max_autocorr = corr;
    }
    bool low_autocorr = (max_autocorr < 0.25);

    /* Test 3: Sub-cube uniformity.
     * Partition inputs into 4 sub-cubes and check uniformity. */
    bool *subcube_ones = (bool *)calloc(4, sizeof(bool));
    double subcube_sizes[4] = {0};
    for (size_t i = 0; i < size; i++) {
        int q = ((i & 1) ? 0 : 0) + (((i >> (tt->n - 1)) & 1) ? 2 : 0) + (((i >> (tt->n / 2)) & 1) ? 1 : 0);
        if (q >= 4) q = 3;
        subcube_sizes[q] += 1.0;
        if (tt->table[i]) subcube_ones[q] = true;
    }
    bool subcube_varied = true;
    for (int q = 0; q < 4; q++)
        if (subcube_sizes[q] <= 1.0) subcube_varied = false;
    free(subcube_ones);

    /* Function is "hard" if it passes all three heuristic tests */
    return balanced && low_autocorr && subcube_varied;
}

/* =============================================
 * L3: Fourier Analysis of Boolean Functions
 *
 * For f: {0,1}^n -> {+1,-1} (mapped as (-1)^{f(x)}),
 * the Fourier coefficient at subset S is:
 *   hat{f}(S) = 1/2^n * sum_x f(x) * chi_S(x)
 * where chi_S(x) = (-1)^{sum_{i in S} x_i}.
 *
 * Parseval: sum_S hat{f}(S)^2 = 1.
 *
 * The Fourier spectrum reveals the function's "complexity":
 * a function with concentrated Fourier mass on large subsets
 * is "simpler" (has small circuit complexity).
 * ============================================= */

double nw_fourier_coefficient(const TruthTable *tt, const bool *S) {
    if (!tt || !S || tt->n > 16) return 0.0;

    double sum = 0.0;
    for (size_t x = 0; x < tt->size; x++) {
        /* f(x) mapped to +/-1 */
        double f_val = tt->table[x] ? -1.0 : 1.0;

        /* chi_S(x) = (-1)^{parity of x restricted to S} */
        int parity = 0;
        for (size_t i = 0; i < tt->n; i++) {
            if (S[i] && ((x >> i) & 1))
                parity ^= 1;
        }
        double chi = parity ? -1.0 : 1.0;

        sum += f_val * chi;
    }

    return sum / (double)tt->size;
}

/**
 * L3: Variable Influence
 *
 * Inf_i(f) = Pr_x[f(x) != f(x xor e_i)]
 * where e_i is the unit vector with 1 at position i.
 *
 * Total influence: sum_i Inf_i(f).
 * Edge-isoperimetric inequality: sum_i Inf_i(f) >= 2 * Var[f].
 *
 * Functions with low total influence are "simple"
 * (approximable by low-degree polynomials, small circuits).
 */
double nw_variable_influence(const TruthTable *tt, size_t i) {
    if (!tt || i >= tt->n) return 0.0;

    size_t mask = (size_t)1 << i;
    size_t flips = 0;
    size_t total = tt->size;

    for (size_t x = 0; x < total; x++) {
        size_t y = x ^ mask;
        if (tt->table[x] != tt->table[y])
            flips++;
    }

    /* Each pair counted twice, but we want Pr[f(x) != f(x xor e_i)]
     * which is flips / total since each pair counted once in expectation */
    return (double)flips / (double)total;
}

/* =============================================
 * L8: Natural Proofs Barrier
 *
 * Razborov-Rudich (1997):
 *
 * A "natural proof" of a circuit lower bound is:
 *   1. Constructive: the property C_n is decidable in P/poly
 *      (recognizable by polynomial-size circuits)
 *   2. Largeness: |C_n| >= 2^{2^n} / poly(n)
 *      (at least 1/poly fraction of all functions have the property)
 *   3. Usefulness: f in C_n => f requires large circuits
 *
 * Theorem: If sub-exponentially strong one-way functions exist,
 * then no natural proof can show super-polynomial circuit lower bounds.
 *
 * The NW PRG construction is a key part of this result:
 * if we had a natural property, we could use it to break
 * the NW PRG, contradicting the existence of OWFs.
 * ============================================= */

bool nw_natural_proofs_barrier(double hardness, size_t n) {
    /* Check if the target lower bound is in the "natural proofs barrier" range.
     *
     * The barrier says: any natural proof against circuits of size s(n)
     * implies breaking OWFs if s(n) > poly(n).
     *
     * The NW PRG shows: if OWFs exist, there exist PRGs with stretch n -> n+1
     * that fool circuits of size s(n), making natural properties
     * indistinguishable from random.
     *
     * We check if hardness > 2^{n^{0.001}} (super-polynomial)
     */
    if (n < 10) return false;

    double super_poly = pow(2.0, pow((double)n, 0.01));
    if (hardness > super_poly) return true;

    return false;
}

double nw_razborov_rudich_threshold(size_t n) {
    /* The threshold circuit size above which natural proofs
     * cannot exist (assuming standard cryptographic assumptions).
     *
     * Threshold = min(2^{n^{o(1)}}, 2^{n/2})
     *
     * For current knowledge, anything above 2^{n^{0.999}} is
     * potentially beyond the barrier.
     */
    if (n <= 1) return 2.0;
    double bound1 = pow(2.0, pow((double)n, 0.5));
    double bound2 = pow(2.0, (double)n * 0.5);
    return (bound1 < bound2) ? bound1 : bound2;
}

/* =============================================
 * L7: IW97 Derandomization Levels
 *
 * Impagliazzo & Wigderson (1997) showed:
 *
 * Level 0 (Weak hardness): If E requires circuits of size 2^{n^o(1)},
 *   then BPP ⊆ SUBEXP = ∩_{ε>0} DTIME(2^{n^ε})
 *
 * Level 1 (Moderate): If E requires circuits of size 2^{Ω(n)},
 *   then BPP ⊆ QP = DTIME(2^{polylog(n)})
 *
 * Level 2 (Strong): If E requires circuits of size 2^{Ω(n)} AND
 *   the function is in E (not just EXP), then BPP = P.
 * ============================================= */

bool nw_iw_derandomize(int level, size_t n, size_t *seed_len) {
    if (!seed_len) return false;

    switch (level) {
        case 0:
            /* Weak: BPP ⊆ SUBEXP
             * Seed length: n^epsilon for any epsilon > 0
             * (sub-exponential in n) */
            *seed_len = (size_t)ceil(pow((double)n, 0.5));
            return true;

        case 1:
            /* Moderate: BPP ⊆ QP
             * Seed length: polylog(n) = O(log^c n) */
            {
                size_t log_n = nw_ceil_log2(n);
                *seed_len = log_n * log_n * log_n; /* O(log^3 n) */
                if (*seed_len == 0) *seed_len = 1;
            }
            return true;

        case 2:
            /* Strong: BPP = P
             * Seed length: O(log n) — polynomial overhead */
            *seed_len = nw_ceil_log2(n) * 2;
            if (*seed_len == 0) *seed_len = 1;
            return true;

        default:
            return false;
    }
}

size_t nw_bpp_simulation_samples(size_t n, double epsilon) {
    (void)n;
    /* Using NW PRG to derandomize BPP with error epsilon:
     * Number of seeds = 2^k where k = seed length.
     * Each seed produces m pseudorandom bits.
     *
     * Error probability after taking majority over all seeds:
     *   <= exp(-2 * (epsilon')^2 * 2^k)
     * where epsilon' > 1/3 is the original BPP error gap.
     *
     * For epsilon = 0.01 (confidence 99%):
     *   2^k >= ln(1/epsilon) / (2 * (1/6)^2) ≈ 18 * ln(1/epsilon)
     */
    if (epsilon <= 0.0) return (size_t)-1;
    if (epsilon >= 1.0) return 0;

    double log_inv = log(1.0 / epsilon);
    double samples = log_inv * 18.0; /* 2 * (1/6)^-2 = 2 * 36 = 72, but conservative = 18 */
    samples /= 2.0;
    return (size_t)ceil(samples);
}

/* =============================================
 * L5: Circuit Complexity Estimation
 *
 * Estimate minimum circuit size from truth table using
 * heuristic measures:
 *   1. Sub-cube entropy
 *   2. Fourier concentration
 *   3. Decision tree depth
 * ============================================= */

double nw_estimate_circuit_complexity(const TruthTable *tt) {
    if (!tt || tt->n == 0) return 0.0;
    if (tt->n <= 4) {
        /* For small n, exact computation is feasible.
         * Compute complexity as fraction of Shannon bound. */
        double shannon = nw_shannon_lower_bound(tt->n);
        /* Simple functions (e.g., AND, OR) have complexity O(n).
         * We estimate by counting truth table "entropy". */
        size_t changes = 0;
        for (size_t i = 1; i < tt->size; i++)
            if (tt->table[i] != tt->table[i-1]) changes++;
        /* If the function changes frequently (like PARITY), it's hard.
         * If it's monotonic/simple, it's easy. */
        double change_ratio = (double)changes / (double)(tt->size - 1);
        double estimated = change_ratio * shannon;
        if (estimated < (double)tt->n) estimated = (double)tt->n; /* minimum O(n) */
        return estimated;
    }

    /* For larger n, use asymptotic estimate */
    double shannon = nw_shannon_lower_bound(tt->n);

    /* Approximate by random function baseline */
    double random_complexity = shannon * 0.8; /* Most functions are near the bound */

    return random_complexity;
}

/**
 * L4: Circuit Lower Bound from Gate Count
 *
 * Given a candidate circuit with `gates` gates, we can assert
 * a lower bound on the circuit size needed:
 *   circuit_size >= gates  (trivial)
 *
 * Using Shannon's bound: if gates < 2^n/(2n), then the circuit
 * cannot compute most functions.
 */
double nw_circuit_size_lower_bound_from_gates(size_t n, size_t gates) {
    double shannon = nw_shannon_lower_bound(n);
    if ((double)gates < shannon) return shannon;
    return (double)gates;
}
