#include "next_bit_tests.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* =========================================================================
 * L1/L5: NIST Statistical Test Suite for Random Number Generators
 *
 * The NIST SP 800-22 test suite evaluates whether a binary sequence
 * is statistically random. These tests are fundamental to assessing
 * whether a PRG's output is computationally indistinguishable from
 * truly random bits.
 *
 * Each test computes a P-value. If P-value >= alpha (typically 0.01),
 * the sequence passes the test; otherwise it is considered non-random.
 *
 * Reference: NIST SP 800-22 Rev 1a (2010)
 * Course: MIT 6.857, Stanford CS255 S8, CMU 15-855 S6
 * ========================================================================= */

/**
 * nist_frequency_monobit - Frequency (Monobit) Test
 *
 * Tests whether the proportion of 1s is approximately 1/2.
 * Test statistic: s_obs = |S_n| / sqrt(n) where S_n = sum(2*bit-1)
 * P-value = erfc(s_obs / sqrt(2))
 *
 * Minimum recommended sequence length: 100 bits.
 *
 * @param bits  Sequence of bits
 * @param n     Number of bits
 * @return      P-value in [0,1]; higher = more random
 *
 * Complexity: O(n) time
 */
double nist_frequency_monobit(const bool *bits, size_t n) {
    if (bits == NULL || n == 0) return 0.0;
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum += bits[i] ? 1.0 : -1.0;
    }
    double s_obs = fabs(sum) / sqrt((double)n);
    double p_value = erfc(s_obs / sqrt(2.0));
    return p_value;
}

/**
 * nist_frequency_block - Frequency Test within a Block
 *
 * Tests whether the proportion of 1s within M-bit blocks is approximately M/2.
 * Chi-squared test on the distribution of 1s counts in N = floor(n/M) blocks.
 *
 * Test statistic: chi^2_obs = 4*M * sum_{i=1..N} (pi_i - 1/2)^2
 * where pi_i = proportion of 1s in block i.
 * P-value = igamc(N/2, chi^2_obs/2) (incomplete gamma function)
 *
 * @param bits       Sequence of bits
 * @param n          Total sequence length
 * @param block_size Size of each block M
 * @return           P-value
 */
double nist_frequency_block(const bool *bits, size_t n, size_t block_size) {
    if (bits == NULL || n == 0 || block_size == 0) return 0.0;
    size_t N = n / block_size;
    if (N == 0) return 0.0;

    double chi_sq = 0.0;
    for (size_t i = 0; i < N; i++) {
        size_t ones = 0;
        for (size_t j = 0; j < block_size; j++) {
            if (bits[i * block_size + j]) ones++;
        }
        double pi = (double)ones / (double)block_size;
        double diff = pi - 0.5;
        chi_sq += diff * diff;
    }
    chi_sq *= 4.0 * (double)block_size;

    /* Approximate incomplete gamma for chi-squared */
    /* P-value = Q(N/2, chi_sq/2) where Q is regularized gamma */
    /* For large N, chi-squared distribution: P-value = 1 - CDF */
    double x = chi_sq / 2.0;
    double a = (double)N / 2.0;

    /* Simple approximation for chi-squared P-value */
    /* Use Wilson-Hilferty transformation for moderate N */
    if (N >= 5) {
        double z = (pow(x / a, 1.0/3.0) - (1.0 - 2.0/(9.0*a))) / sqrt(2.0/(9.0*a));
        double p = 0.5 * erfc(z / sqrt(2.0));
        return p;
    }

    /* For small N, use sum of Poisson probabilities */
    /* Approximation: P-value = 1 - sum_{k=0}^{N/2-1} (x^k/k!) * e^{-x} */
    if (x < 0.0) return 0.0;
    double sum = 0.0, term = 1.0;
    for (size_t k = 0; k < (size_t)a && term > 1e-15; k++) {
        sum += term;
        term *= x / (double)(k + 1);
    }
    double p_value = 1.0 - sum * exp(-x);
    if (p_value < 0.0) p_value = 0.0;
    if (p_value > 1.0) p_value = 1.0;
    return p_value;
}

/**
 * nist_runs_test - Runs Test
 *
 * Tests whether the number of runs (uninterrupted sequences of identical bits)
 * is consistent with randomness.
 *
 * A run is a maximal consecutive subsequence of 1s or 0s.
 * For random sequence of length n with pi proportion of 1s:
 *   E[runs] = 2*n*pi*(1-pi) + 1
 *   Var[runs] = 2*n*pi*(1-pi)*(2*n*pi*(1-pi)-1) / (n-1)
 *
 * @param bits  Sequence of bits
 * @param n     Number of bits
 * @return      P-value
 */
double nist_runs_test(const bool *bits, size_t n) {
    if (bits == NULL || n < 2) return 0.0;

    /* Count proportion of 1s */
    size_t ones = 0;
    for (size_t i = 0; i < n; i++) if (bits[i]) ones++;
    double pi = (double)ones / (double)n;

    /* Frequency pre-test: tau = 2/sqrt(n) */
    double tau = 2.0 / sqrt((double)n);
    if (fabs(pi - 0.5) >= tau) return 0.0;

    /* Count runs */
    size_t runs = 1;
    for (size_t i = 1; i < n; i++) {
        if (bits[i] != bits[i-1]) runs++;
    }

    double expected = 2.0 * (double)n * pi * (1.0 - pi) + 1.0;
    double variance = 2.0 * (double)n * pi * (1.0 - pi) *
                      (2.0 * (double)n * pi * (1.0 - pi) - 1.0) /
                      ((double)n - 1.0);

    if (variance <= 0.0) return 0.0;
    double z = ((double)runs - expected) / sqrt(variance);
    double p_value = erfc(fabs(z) / sqrt(2.0));
    return p_value;
}

/**
 * nist_longest_run_ones - Longest Run of Ones Test
 *
 * Tests whether the longest run of 1s within M-bit blocks is consistent
 * with randomness. Based on the distribution of longest runs in random
 * binary sequences.
 *
 * @param bits  Sequence of bits
 * @param n     Number of bits (must be >= 128)
 * @return      P-value
 */
double nist_longest_run_ones(const bool *bits, size_t n) {
    if (bits == NULL || n < 128) return 0.0;

    /* For n >= 128: M=8, N=16 blocks, K=3 categories */
    size_t M = 8, N = n / 8;
    if (n >= 6272) { M = 128; N = n / 128; }
    else if (n >= 750000) { M = 10000; N = n / 10000; }

    /* For M=8: categories: <=1, 2, 3, >=4 */
    /* Expected probabilities for longest run in 8-bit block */
    double pi_arr[] = {0.2148, 0.3672, 0.2305, 0.1875};
    size_t K = 4;
    size_t counts[4] = {0};

    for (size_t b = 0; b < N; b++) {
        size_t start = b * M;
        size_t longest = 0, current = 0;
        for (size_t i = 0; i < M; i++) {
            if (bits[start + i]) {
                current++;
                if (current > longest) longest = current;
            } else {
                current = 0;
            }
        }
        if (longest <= 1) counts[0]++;
        else if (longest == 2) counts[1]++;
        else if (longest == 3) counts[2]++;
        else counts[3]++;
    }

    /* Chi-squared test */
    double chi_sq = 0.0;
    for (size_t i = 0; i < K; i++) {
        double expected = pi_arr[i] * (double)N;
        double diff = (double)counts[i] - expected;
        chi_sq += diff * diff / expected;
    }

    /* Approximate P-value from chi-squared with K-1 degrees of freedom */
    double p_value = exp(-chi_sq / 2.0);
    if (p_value > 1.0) p_value = 1.0;
    return p_value;
}

/**
 * nist_serial_test - Serial Test (dependence test)
 *
 * Tests the uniformity of the distribution of all possible overlapping
 * m-bit patterns. A random sequence should have each m-bit pattern
 * appear with roughly equal frequency.
 *
 * Uses the psi-squared statistic based on transition counts.
 *
 * @param bits  Sequence of bits
 * @param n     Sequence length
 * @param m     Pattern length (typically <= floor(log2(n)) - 2)
 * @return      P-value
 */
double nist_serial_test(const bool *bits, size_t n, size_t m) {
    if (bits == NULL || n == 0 || m == 0 || m > 16) return 0.0;

    /* Count overlapping m-bit, (m-1)-bit, and (m-2)-bit patterns */
    size_t patterns_m = (size_t)1 << m;
    size_t patterns_m1 = (size_t)1 << (m - 1);
    size_t patterns_m2 = (m >= 2) ? (size_t)1 << (m - 2) : 1;

    size_t *counts_m = (size_t *)calloc(patterns_m, sizeof(size_t));
    size_t *counts_m1 = (size_t *)calloc(patterns_m1, sizeof(size_t));
    size_t *counts_m2 = (size_t *)calloc(patterns_m2, sizeof(size_t));

    if (counts_m == NULL || counts_m1 == NULL || counts_m2 == NULL) {
        free(counts_m); free(counts_m1); free(counts_m2);
        return 0.0;
    }

    /* Augment sequence: append first m-1 bits */
    size_t aug_n = n + m - 1;
    bool *aug = (bool *)malloc(aug_n * sizeof(bool));
    if (aug == NULL) { free(counts_m); free(counts_m1); free(counts_m2); return 0.0; }
    memcpy(aug, bits, n * sizeof(bool));
    memcpy(aug + n, bits, (m-1) * sizeof(bool));

    /* Count m-bit patterns */
    for (size_t i = 0; i < n; i++) {
        size_t idx = 0;
        for (size_t j = 0; j < m; j++) {
            idx = (idx << 1) | (aug[i+j] ? 1 : 0);
        }
        if (idx < patterns_m) counts_m[idx]++;
    }

    /* Count (m-1)-bit patterns */
    for (size_t i = 0; i < n; i++) {
        size_t idx = 0;
        for (size_t j = 0; j < m - 1; j++) {
            idx = (idx << 1) | (aug[i+j] ? 1 : 0);
        }
        if (idx < patterns_m1) counts_m1[idx]++;
    }

    /* Count (m-2)-bit patterns */
    if (m >= 2) {
        for (size_t i = 0; i < n; i++) {
            size_t idx = 0;
            for (size_t j = 0; j < m - 2; j++) {
                idx = (idx << 1) | (aug[i+j] ? 1 : 0);
            }
            if (idx < patterns_m2) counts_m2[idx]++;
        }
    }

    /* psi^2_m = (2^m / n) * sum(counts^2) - n */
    double psi_sq_m = 0.0, psi_sq_m1 = 0.0, psi_sq_m2 = 0.0;

    for (size_t i = 0; i < patterns_m; i++)
        psi_sq_m += (double)counts_m[i] * (double)counts_m[i];
    psi_sq_m = psi_sq_m * (double)patterns_m / (double)n - (double)n;

    for (size_t i = 0; i < patterns_m1; i++)
        psi_sq_m1 += (double)counts_m1[i] * (double)counts_m1[i];
    psi_sq_m1 = psi_sq_m1 * (double)patterns_m1 / (double)n - (double)n;

    if (m >= 2) {
        for (size_t i = 0; i < patterns_m2; i++)
            psi_sq_m2 += (double)counts_m2[i] * (double)counts_m2[i];
        psi_sq_m2 = psi_sq_m2 * (double)patterns_m2 / (double)n - (double)n;
    }

    /* delta_psi^2_m = psi^2_m - psi^2_{m-1} */
    double delta1 = psi_sq_m - psi_sq_m1;
    double delta2 = (m >= 2) ? psi_sq_m1 - psi_sq_m2 : psi_sq_m1;

    /* P-value from chi-squared */
    size_t df = (size_t)1 << (m - 1); /* degrees of freedom */
    double p_value1 = (delta1 > 0) ? exp(-delta1 / 2.0) : 1.0;
    double p_value2 = (delta2 > 0) ? exp(-delta2 / 2.0) : 1.0;

    free(counts_m); free(counts_m1); free(counts_m2); free(aug);
    return (p_value1 < p_value2) ? p_value1 : p_value2;
}

/**
 * nist_approximate_entropy - Approximate Entropy Test
 *
 * Compares the frequency of overlapping m-bit and (m+1)-bit patterns.
 * The approximate entropy ApEn(m) = Phi^m - Phi^{m+1} where
 * Phi^m = sum_i (count_i/n) * log(count_i/n).
 *
 * For random sequences, ApEn should be close to ln(2) for large m.
 *
 * @param bits  Sequence
 * @param n     Length
 * @param m     Block size (typically m < floor(log2(n)) - 5)
 * @return      P-value
 */
double nist_approximate_entropy(const bool *bits, size_t n, size_t m) {
    if (bits == NULL || n == 0 || m == 0 || m > 8) return 0.0;

    size_t patterns = (size_t)1 << (m + 1);
    size_t *counts = (size_t *)calloc(patterns, sizeof(size_t));
    if (counts == NULL) return 0.0;

    /* Count overlapping (m+1)-bit patterns */
    bool *aug = (bool *)malloc((n + m) * sizeof(bool));
    if (aug == NULL) { free(counts); return 0.0; }
    memcpy(aug, bits, n * sizeof(bool));
    memcpy(aug + n, bits, m * sizeof(bool));

    size_t mask = patterns - 1;
    size_t current = 0;
    /* Initialize with first m bits */
    for (size_t i = 0; i < m; i++) {
        current = ((current << 1) | (aug[i] ? 1 : 0)) & mask;
    }

    for (size_t i = 0; i < n; i++) {
        current = ((current << 1) | (aug[i + m] ? 1 : 0)) & mask;
        counts[current]++;
    }

    /* Compute Phi^{m+1} */
    double phi = 0.0;
    for (size_t i = 0; i < patterns; i++) {
        if (counts[i] > 0) {
            double p = (double)counts[i] / (double)n;
            phi += p * log(p);
        }
    }

    /* Approximate P-value: chi-squared with 2^m degrees of freedom */
    double chi_sq = 2.0 * (double)n * (log(2.0) + phi);
    double p_value = exp(-chi_sq / 4.0); /* simplified */
    if (p_value > 1.0) p_value = 1.0;
    if (p_value < 0.0) p_value = 0.0;

    free(counts); free(aug);
    return p_value;
}

/**
 * nist_cumulative_sums - Cumulative Sums (Cusum) Test
 *
 * Tests whether the cumulative sum of adjusted bits (-1,+1) stays
 * near zero. Maximum excursion from zero is the test statistic.
 *
 * @param bits    Sequence
 * @param n       Length
 * @param forward true = forward direction, false = reverse
 * @return        P-value
 */
double nist_cumulative_sums(const bool *bits, size_t n, bool forward) {
    if (bits == NULL || n == 0) return 0.0;

    double max_z = 0.0;
    double sum = 0.0;

    for (size_t i = 0; i < n; i++) {
        size_t idx = forward ? i : (n - 1 - i);
        sum += bits[idx] ? 1.0 : -1.0;
        double abs_sum = (sum < 0.0) ? -sum : sum;
        if (abs_sum > max_z) max_z = abs_sum;
    }

    /* P-value approximation for the maximum partial sum */
    double z = max_z / sqrt((double)n);
    double p_value = 0.0;
    /* Summation over k = (-n/z - 1)/2 to (n/z - 1)/2 */
    double term1 = erfc(z / sqrt(2.0));
    /* Simplified: the standard normal tail approximation */
    p_value = term1;

    if (p_value > 1.0) p_value = 1.0;
    return p_value;
}

/* =========================================================================
 * Additional Statistical Tests
 * ========================================================================= */

/**
 * chi_squared_test - Pearson's Chi-Squared Test
 *
 * Tests goodness-of-fit between observed and expected frequencies.
 * chi^2 = sum (O_i - E_i)^2 / E_i
 *
 * @param observed  Observed counts per bucket
 * @param expected  Expected counts per bucket
 * @param buckets   Number of buckets
 * @return          P-value (approximate)
 */
double chi_squared_test(const size_t *observed, const double *expected, size_t buckets) {
    if (observed == NULL || expected == NULL || buckets == 0) return 0.0;

    double chi_sq = 0.0;
    for (size_t i = 0; i < buckets; i++) {
        if (expected[i] > 0.0) {
            double diff = (double)observed[i] - expected[i];
            chi_sq += diff * diff / expected[i];
        }
    }

    /* Approximate P-value for chi-squared distribution */
    double df = (double)(buckets - 1);
    if (df <= 0.0) df = 1.0;
    double p_value = exp(-chi_sq / (2.0 * sqrt(df)));
    if (p_value > 1.0) p_value = 1.0;
    return p_value;
}

/**
 * kolmogorov_smirnov_test - Kolmogorov-Smirnov Test
 *
 * Tests whether a sample comes from a uniform distribution.
 * D_n = max_i |F_n(x_i) - F(x_i)| where F_n is the empirical CDF.
 *
 * @param samples  Array of sample values in [0,1]
 * @param n        Number of samples
 * @return         P-value (Kolmogorov approximation)
 */
double kolmogorov_smirnov_test(const double *samples, size_t n) {
    if (samples == NULL || n == 0) return 0.0;

    /* Sort samples (insertion sort for simplicity) */
    double *sorted = (double *)malloc(n * sizeof(double));
    if (sorted == NULL) return 0.0;
    memcpy(sorted, samples, n * sizeof(double));

    for (size_t i = 1; i < n; i++) {
        double key = sorted[i];
        size_t j = i;
        while (j > 0 && sorted[j-1] > key) {
            sorted[j] = sorted[j-1];
            j--;
        }
        sorted[j] = key;
    }

    /* Compute KS statistic D_n */
    double d_max = 0.0;
    for (size_t i = 0; i < n; i++) {
        /* F(x) = x for uniform[0,1]; F_n(x_i) = (i+1)/n */
        double fn = (double)(i + 1) / (double)n;
        double d1 = fn - sorted[i];
        double d2 = sorted[i] - (double)i / (double)n;
        if (d1 < 0.0) d1 = -d1;
        if (d2 < 0.0) d2 = -d2;
        if (d1 > d_max) d_max = d1;
        if (d2 > d_max) d_max = d2;
    }

    free(sorted);

    /* Kolmogorov approximation: P(D_n > lambda) ~ 2 * sum_{k=1}^{inf} (-1)^{k-1} e^{-2k^2 lambda^2 */
    double lambda = (sqrt((double)n) + 0.12 + 0.11 / sqrt((double)n)) * d_max;
    double p_value = 2.0 * exp(-2.0 * lambda * lambda);
    if (p_value > 1.0) p_value = 1.0;
    return p_value;
}

/**
 * binary_entropy - Shannon Entropy of binary sequence
 *
 * H = -p * log2(p) - (1-p) * log2(1-p) where p = proportion of 1s.
 * Maximum entropy = 1 bit/bit at p=0.5.
 *
 * @param bits  Sequence
 * @param n     Length
 * @return      Entropy in bits (0 to 1)
 */
double binary_entropy(const bool *bits, size_t n) {
    if (bits == NULL || n == 0) return 0.0;
    size_t ones = 0;
    for (size_t i = 0; i < n; i++) if (bits[i]) ones++;
    double p = (double)ones / (double)n;
    if (p <= 0.0 || p >= 1.0) return 0.0;
    return -p * log2(p) - (1.0 - p) * log2(1.0 - p);
}

/**
 * linear_complexity - Linear complexity via Berlekamp-Massey algorithm
 *
 * Linear complexity = length of shortest LFSR that generates the sequence.
 * For a truly random sequence of length n, expected linear complexity ~ n/2.
 *
 * Uses the Berlekamp-Massey algorithm over GF(2).
 * Complexity: O(n^2) time, O(n) space.
 *
 * @param bits  Sequence
 * @param n     Length
 * @return      Linear complexity (0 to n)
 */
size_t linear_complexity(const bool *bits, size_t n) {
    if (bits == NULL || n == 0) return 0;

    size_t *c = (size_t *)calloc(n + 1, sizeof(size_t));
    size_t *b = (size_t *)calloc(n + 1, sizeof(size_t));
    if (c == NULL || b == NULL) { free(c); free(b); return 0; }

    c[0] = 1; b[0] = 1;
    size_t L = 0, m = -1;
    size_t mm = 0;

    for (size_t i = 0; i < n; i++) {
        /* Compute discrepancy d */
        size_t d = bits[i] ? 1 : 0;
        for (size_t j = 1; j <= L && j <= i; j++) {
            if (c[j] && bits[i - j]) d ^= 1;
        }

        if (d == 1) {
            /* Save copy of c */
            size_t *t = (size_t *)calloc(n + 1, sizeof(size_t));
            if (t == NULL) { free(c); free(b); return L; }
            memcpy(t, c, (n+1) * sizeof(size_t));

            /* c = c - d * d^{-1} * x^{i-m} * b */
            /* Over GF(2): d = d^{-1} = 1, so c = c + x^{i-m} * b */
            size_t shift = i - mm;
            for (size_t j = 0; j <= n; j++) {
                if (j + shift <= n && b[j]) {
                    c[j + shift] ^= 1;
                }
            }

            if (2 * L <= i) {
                L = i + 1 - L;
                mm = i;
                memcpy(b, t, (n+1) * sizeof(size_t));
            }
            free(t);
        }
    }

    free(c); free(b);
    return L;
}

/**
 * autocorrelation - Autocorrelation at given lag
 *
 * R(lag) = (1/(n-lag)) * sum_{i=0}^{n-lag-1} (2b_i-1)(2b_{i+lag}-1)
 *
 * For random bits: E[R(lag)] ~ 0, Var[R(lag)] ~ 1/n.
 *
 * @param bits  Sequence
 * @param n     Length
 * @param lag   Lag distance
 * @return      Autocorrelation coefficient in [-1,1]
 */
double autocorrelation(const bool *bits, size_t n, size_t lag) {
    if (bits == NULL || n == 0 || lag >= n) return 0.0;

    double sum = 0.0;
    size_t count = n - lag;
    for (size_t i = 0; i < count; i++) {
        double a = bits[i] ? 1.0 : -1.0;
        double b = bits[i + lag] ? 1.0 : -1.0;
        sum += a * b;
    }
    return sum / (double)count;
}

/* =========================================================================
 * Combined test battery
 * ========================================================================= */

/**
 * next_bit_test_battery - Run a battery of next-bit predictability tests.
 *
 * Computes multiple statistical tests and returns a composite score.
 * The composite is the minimum P-value across all tests.
 *
 * @param bits    Sequence to test
 * @param n       Sequence length
 * @param results Array of NISTTestResult to fill
 * @param max_tests Max results to fill
 * @return        Minimum P-value across tests
 */
double next_bit_test_battery(const bool *bits, size_t n,
                              NISTTestResult *results, size_t max_tests) {
    if (bits == NULL || results == NULL || n == 0 || max_tests == 0) return 0.0;

    size_t idx = 0;
    double min_p = 1.0;

    #define RUN_TEST(ttype, p_val) do { \
        if (idx < max_tests) { \
            results[idx].type = ttype; \
            results[idx].p_value = p_val; \
            results[idx].passed = (p_val >= 0.01); \
            idx++; \
        } \
        if (p_val < min_p) min_p = p_val; \
    } while(0)

    RUN_TEST(TEST_FREQUENCY_MONOBIT, nist_frequency_monobit(bits, n));
    if (n >= 100) {
        RUN_TEST(TEST_FREQUENCY_BLOCK, nist_frequency_block(bits, n, n/10 > 0 ? n/10 : 1));
    }
    RUN_TEST(TEST_RUNS, nist_runs_test(bits, n));
    if (n >= 128) {
        RUN_TEST(TEST_LONGEST_RUN, nist_longest_run_ones(bits, n));
    }
    if (n >= 64) {
        size_t m = 2;
        while ((size_t)1 << m <= n) m++;
        m = (m > 2) ? m - 2 : 2;
        RUN_TEST(TEST_SERIAL, nist_serial_test(bits, n, m));
    }
    if (n >= 64) {
        size_t m = 2;
        while ((size_t)1 << (m+1) <= n) m++;
        m = (m > 1) ? m - 1 : 1;
        if (m <= 8) {
            RUN_TEST(TEST_APPROX_ENTROPY, nist_approximate_entropy(bits, n, m));
        }
    }
    RUN_TEST(TEST_CUMULATIVE_SUMS, nist_cumulative_sums(bits, n, true));

    #undef RUN_TEST

    return min_p;
}
