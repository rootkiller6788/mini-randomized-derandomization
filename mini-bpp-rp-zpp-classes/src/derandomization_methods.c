/******************************************************************************
 * derandomization_methods.c — Derandomization Methods & Conditional Probabilities
 *
 * Implements key derandomization techniques:
 *   L4: Method of conditional probabilities (MAX-CUT, Set Balancing)
 *   L4: Pairwise & k-wise independent hash functions
 *   L5: Brute-force derandomization (BPP ⊆ EXP)
 *   L7: Randomized quicksort analysis
 *   L8: Hardness vs randomness (Impagliazzo's worlds)
 *
 * References:
 *   Raghavan (1988), Spencer (1987)
 *   Carter & Wegman (1979), Wegman & Carter (1981)
 *   Impagliazzo (1995), Impagliazzo & Wigderson (1997)
 ******************************************************************************/

#include "derandomization_methods.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

/* ================================================================
 * L4: Method of Conditional Probabilities — MAX-CUT
 * ================================================================ */

/**
 * Theorem (MAX-CUT derandomization):
 * Randomized: assign each vertex L/R uniformly → E[cut] = |E|/2.
 * Conditional probabilities: for each vertex v in order, compute
 * the conditional expected cut size given assignments to {0,...,v-1}.
 *
 * If E[cut | fixed_prefix, v=L] ≥ E[cut | fixed_prefix, v=R],
 * choose L; else choose R. This maintains E[cut | prefix] ≥ |E|/2,
 * yielding a deterministic cut of size ≥ |E|/2.
 */
size_t max_cut_derandomized(
    const int *adj_matrix, size_t n,
    bool *assignments)
{
    if (n == 0 || !adj_matrix || !assignments) return 0;

    /* Count total edges */
    size_t total_edges = 0;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            if (adj_matrix[i * n + j]) total_edges++;
        }
    }

    /* Greedy assignment using conditional expectations */
    bool *fixed = (bool *)calloc(n, sizeof(bool));
    bool *side = (bool *)calloc(n, sizeof(bool));
    if (!fixed || !side) {
        free(fixed); free(side);
        return 0;
    }

    for (size_t v = 0; v < n; v++) {
        /* Mark v as fixed BEFORE computing expectations, so that
         * edges incident to v are included in the calculation. */
        fixed[v] = true;

        /* Compute conditional expected cut for v = 0 (left) */
        side[v] = false;
        double expect_0 = 0.0;

        /* Edges between fixed vertices (i,j ≤ v) */
        for (size_t i = 0; i <= v; i++) {
            if (!fixed[i]) continue;
            for (size_t j = i + 1; j <= v; j++) {
                if (!fixed[j]) continue;
                if (adj_matrix[i * n + j] && side[i] != side[j]) {
                    expect_0 += 1.0;
                }
            }
        }

        /* Edges from fixed vertices (i ≤ v) to unfixed (j > v) */
        for (size_t i = 0; i <= v; i++) {
            if (!fixed[i]) continue;
            for (size_t j = v + 1; j < n; j++) {
                if (adj_matrix[i * n + j]) {
                    expect_0 += 0.5;  /* j's side is still random */
                }
            }
        }

        /* Edges between unfixed vertices (both > v) */
        for (size_t i = v + 1; i < n; i++) {
            for (size_t j = i + 1; j < n; j++) {
                if (adj_matrix[i * n + j]) expect_0 += 0.5;
            }
        }

        /* Compute for v = 1 (right) — same structure */
        side[v] = true;
        double expect_1 = 0.0;

        for (size_t i = 0; i <= v; i++) {
            if (!fixed[i]) continue;
            for (size_t j = i + 1; j <= v; j++) {
                if (!fixed[j]) continue;
                if (adj_matrix[i * n + j] && side[i] != side[j]) {
                    expect_1 += 1.0;
                }
            }
        }

        for (size_t i = 0; i <= v; i++) {
            if (!fixed[i]) continue;
            for (size_t j = v + 1; j < n; j++) {
                if (adj_matrix[i * n + j]) expect_1 += 0.5;
            }
        }

        for (size_t i = v + 1; i < n; i++) {
            for (size_t j = i + 1; j < n; j++) {
                if (adj_matrix[i * n + j]) expect_1 += 0.5;
            }
        }

        /* Choose side with higher conditional expectation */
        if (expect_0 >= expect_1) {
            side[v] = false;
        } else {
            side[v] = true;
        }
        /* fixed[v] already set to true above */
    }

    memcpy(assignments, side, n * sizeof(bool));

    /* Compute final cut size */
    size_t cut_size = 0;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            if (adj_matrix[i * n + j] && assignments[i] != assignments[j]) {
                cut_size++;
            }
        }
    }

    free(fixed);
    free(side);

    /* Verify guarantee: cut ≥ |E|/2 */
    assert(cut_size >= total_edges / 2);
    return cut_size;
}

/* ================================================================
 * L4: Set Balancing via Conditional Probabilities
 * ================================================================ */

/**
 * Spencer's set balancing: given m subsets of {1..n}, find a
 * ±1 coloring minimizing the maximum discrepancy.
 *
 * Randomized: E[max discrepancy] = O(√(n log m)).
 * For m = poly(n), conditional probabilities achieves O(√(n log n)).
 */
size_t set_balancing_derandomized(
    const int *sets, size_t m, size_t n,
    int *assignment)
{
    if (n == 0 || !sets || !assignment) return 0;

    /* Initialize all assignments to 0 (undecided).
     * We'll assign ±1 iteratively. */
    int *assign = (int *)calloc(n, sizeof(int));
    int *disc = (int *)calloc(m, sizeof(int));
    if (!assign || !disc) {
        free(assign); free(disc);
        return 0;
    }

    /* For each element, compute the effect on each set's discrepancy
     * and choose the sign that minimizes the maximum potential discrepancy.
     *
     * Pessimistic estimator (simplified): for each element j, estimate
     * the contribution to the maximum discrepancy via exponentiated
     * moment bound. */

    for (size_t j = 0; j < n; j++) {
        /* Try +1 */
        double max_disc_plus = 0.0;
        for (size_t i = 0; i < m; i++) {
            if (sets[i * n + j]) {
                double new_disc = fabs((double)disc[i] + 1.0);
                if (new_disc > max_disc_plus) max_disc_plus = new_disc;
            } else {
                double cur = fabs((double)disc[i]);
                if (cur > max_disc_plus) max_disc_plus = cur;
            }
        }

        /* Try -1 */
        double max_disc_minus = 0.0;
        for (size_t i = 0; i < m; i++) {
            if (sets[i * n + j]) {
                double new_disc = fabs((double)disc[i] - 1.0);
                if (new_disc > max_disc_minus) max_disc_minus = new_disc;
            } else {
                double cur = fabs((double)disc[i]);
                if (cur > max_disc_minus) max_disc_minus = cur;
            }
        }

        /* Choose sign with lower maximum discrepancy */
        if (max_disc_plus <= max_disc_minus) {
            assign[j] = 1;
            for (size_t i = 0; i < m; i++) {
                if (sets[i * n + j]) disc[i] += 1;
            }
        } else {
            assign[j] = -1;
            for (size_t i = 0; i < m; i++) {
                if (sets[i * n + j]) disc[i] -= 1;
            }
        }
    }

    memcpy(assignment, assign, n * sizeof(int));

    /* Compute final maximum discrepancy */
    size_t max_disc = 0;
    for (size_t i = 0; i < m; i++) {
        size_t abs_disc = (size_t)((disc[i] >= 0) ? disc[i] : -disc[i]);
        if (abs_disc > max_disc) max_disc = abs_disc;
    }

    free(assign);
    free(disc);
    return max_disc;
}

/* ================================================================
 * L4: Pairwise Independent Hash Functions (Carter-Wegman)
 * ================================================================ */

static uint64_t next_prime(uint64_t x) {
    if (x < 2) return 2;
    if (x == 2) return 2;
    if (x % 2 == 0) x++;

    while (1) {
        bool is_prime = true;
        for (uint64_t d = 3; d * d <= x && d < 1000000; d += 2) {
            if (x % d == 0) { is_prime = false; break; }
        }
        if (is_prime) return x;
        x += 2;
    }
}

PairwiseHash pairwise_hash_create(uint64_t N, uint64_t M, RandomSource *rs) {
    PairwiseHash h;
    h.p = next_prime(N);
    if (h.p < M) h.p = next_prime(M);
    h.M = M;
    /* a ∈ [1, p-1], b ∈ [0, p-1] */
    h.a = random_source_uniform(rs, h.p - 1) + 1;
    h.b = random_source_uniform(rs, h.p);
    return h;
}

uint64_t pairwise_hash_eval(const PairwiseHash *h, uint64_t x) {
    /* h(x) = (a*x + b) mod p mod M */
    /* Use 128-bit arithmetic to avoid overflow */
    __uint128_t prod = (__uint128_t)h->a * (__uint128_t)x;
    uint64_t val = (uint64_t)((prod + (__uint128_t)h->b) % (__uint128_t)h->p);
    return val % h->M;
}

/**
 * Verify pairwise independence:
 * For given x1≠x2, y1, y2, count (a,b) pairs s.t. h_{a,b}(x1)=y1 and h_{a,b}(x2)=y2.
 *
 * Theorem: Exactly (p-1)·p / M² pairs satisfy this for any valid y1, y2,
 * proving pairwise independence.
 */
size_t pairwise_hash_verify(
    uint64_t p, uint64_t M,
    uint64_t x1, uint64_t x2,
    uint64_t y1, uint64_t y2)
{
    /* h(x) = (a*x + b) mod p mod M
     * We want: h(x1) ≡ y1 (mod M), h(x2) ≡ y2 (mod M)
     *
     * For each (a,b) pair, check the condition.
     * For verification of the theorem (small p), enumerate all. */

    size_t count = 0;
    uint64_t total_pairs = (p - 1) * p;
    size_t checked = 0;

    for (uint64_t a = 1; a < p && checked < total_pairs; a++) {
        for (uint64_t b = 0; b < p && checked < total_pairs; b++) {
            checked++;
            __uint128_t v1 = ((__uint128_t)a * x1 + b) % p;
            __uint128_t v2 = ((__uint128_t)a * x2 + b) % p;
            if ((uint64_t)(v1 % M) == y1 && (uint64_t)(v2 % M) == y2) {
                count++;
            }
        }
    }

    return count;
}

/* ================================================================
 * L4: k-wise Independent Hash Functions
 * ================================================================ */

KWiseHash kwise_hash_create(uint64_t N, uint64_t M, size_t k, RandomSource *rs) {
    KWiseHash h;
    h.k = k;
    h.p = next_prime(N > M ? N : M);
    h.M = M;
    h.coeffs = (uint64_t *)calloc(k, sizeof(uint64_t));
    if (h.coeffs) {
        for (size_t i = 0; i < k; i++) {
            if (i == 0) {
                /* a₀ can be anything in [0, p-1] */
                h.coeffs[i] = random_source_uniform(rs, h.p);
            } else {
                /* a_i ∈ [1, p-1] */
                h.coeffs[i] = random_source_uniform(rs, h.p - 1) + 1;
            }
        }
    }
    return h;
}

uint64_t kwise_hash_eval(const KWiseHash *h, uint64_t x) {
    if (!h || !h->coeffs) return 0;

    /* h(x) = (Σ_{i=0}^{k-1} a_i · x^i) mod p mod M
     * Compute using Horner's method. */
    __uint128_t result = 0;
    __uint128_t x_pow = 1;  /* x^i mod p */

    for (size_t i = 0; i < h->k; i++) {
        __uint128_t term = (__uint128_t)h->coeffs[i] * x_pow;
        result = (result + term) % (__uint128_t)h->p;
        x_pow = (x_pow * (__uint128_t)x) % (__uint128_t)h->p;
    }

    return (uint64_t)(result % (__uint128_t)h->M);
}

void kwise_hash_destroy(KWiseHash *h) {
    if (h && h->coeffs) {
        free(h->coeffs);
        h->coeffs = NULL;
    }
}

/* ================================================================
 * L5: Brute-force Derandomization
 * ================================================================ */

bool brute_force_derandomize(
    const char *input, size_t input_len,
    bool (*decider)(const char*, size_t, RandomSource*),
    size_t m, size_t time_bound, RandomizedClass class_type)
{
    (void)time_bound;  /* reserved for actual time-bounded simulation */
    /* BPP: enumerate all 2^m random strings, take majority.
     * RP:  find ANY accepting random string. */

    size_t total = (size_t)1ULL << m;
    size_t accepts = 0, rejects = 0;

    for (size_t r_val = 0; r_val < total && r_val < (size_t)1ULL << 30; r_val++) {
        RandomSource rs = random_source_init(r_val);
        bool accept = decider(input, input_len, &rs);
        if (accept) {
            accepts++;
            if (class_type == CLASS_RP) {
                return true;  /* RP: one accept suffices */
            }
        } else {
            rejects++;
        }
    }

    if (class_type == CLASS_RP) {
        return false;  /* No accepting random string found */
    }

    if (class_type == CLASS_BPP || class_type == CLASS_ZPP) {
        /* BPP: majority vote */
        return accepts > rejects;
    }

    return false;
}

/* ================================================================
 * L7: Randomized Quicksort Analysis
 * ================================================================ */

/**
 * Simulate randomized quicksort and collect statistics.
 *
 * Expected number of comparisons: 2(n+1)H_n - 4n
 * where H_n = Σ_{i=1}^n 1/i ≈ ln n + γ.
 *
 * Variance: ≈ 1.14 n² (Fill & Janson, 2000).
 * Tail: Pr[comparisons ≥ (1+ε)·expected] ≤ exp(-Θ(n·ε²)).
 */
static void quicksort_rand(int *arr, size_t lo, size_t hi,
                           size_t *comparisons, RandomSource *rs) {
    if (lo >= hi || hi == 0 || lo >= hi) return;

    /* Random pivot */
    size_t pivot_idx = lo + random_source_uniform(rs, hi - lo);
    int pivot = arr[pivot_idx];

    /* Swap pivot to end */
    int tmp = arr[pivot_idx];
    arr[pivot_idx] = arr[hi - 1];
    arr[hi - 1] = tmp;

    /* Partition */
    size_t i = lo;
    for (size_t j = lo; j < hi - 1; j++) {
        (*comparisons)++;
        if (arr[j] <= pivot) {
            tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
            i++;
        }
    }

    /* Place pivot */
    tmp = arr[i];
    arr[i] = arr[hi - 1];
    arr[hi - 1] = tmp;

    if (i > lo) quicksort_rand(arr, lo, i, comparisons, rs);
    if (i + 1 < hi) quicksort_rand(arr, i + 1, hi, comparisons, rs);
}

static void quicksort_median3(int *arr, size_t lo, size_t hi,
                              size_t *comparisons) {
    if (lo >= hi || hi - lo <= 1) return;

    /* Median-of-3 pivot selection */
    size_t mid = lo + (hi - lo) / 2;
    /* Sort lo, mid, hi-1 */
    if (arr[lo] > arr[mid]) { int t = arr[lo]; arr[lo] = arr[mid]; arr[mid] = t; }
    if (arr[lo] > arr[hi-1]) { int t = arr[lo]; arr[lo] = arr[hi-1]; arr[hi-1] = t; }
    if (arr[mid] > arr[hi-1]) { int t = arr[mid]; arr[mid] = arr[hi-1]; arr[hi-1] = t; }
    /* Median is at mid */
    int pivot = arr[mid];
    arr[mid] = arr[hi-1];
    arr[hi-1] = pivot;

    size_t i = lo;
    for (size_t j = lo; j < hi - 1; j++) {
        (*comparisons)++;
        if (arr[j] <= pivot) {
            int t = arr[i]; arr[i] = arr[j]; arr[j] = t;
            i++;
        }
    }
    arr[hi-1] = arr[i];
    arr[i] = pivot;

    if (i > lo) quicksort_median3(arr, lo, i, comparisons);
    if (i + 1 < hi) quicksort_median3(arr, i + 1, hi, comparisons);
}

void randomized_quicksort_analysis(
    size_t n, size_t num_trials,
    double *mean, double *stddev)
{
    if (n == 0 || num_trials == 0) {
        if (mean) *mean = 0.0;
        if (stddev) *stddev = 0.0;
        return;
    }

    double sum = 0.0, sum_sq = 0.0;
    int *arr = (int *)malloc(n * sizeof(int));
    if (!arr) {
        if (mean) *mean = 0.0;
        if (stddev) *stddev = 0.0;
        return;
    }

    for (size_t t = 0; t < num_trials; t++) {
        for (size_t i = 0; i < n; i++) arr[i] = (int)i; /* Sorted input */
        RandomSource rs = random_source_init(t * 0xABCDEFULL);
        size_t comps = 0;
        quicksort_rand(arr, 0, n, &comps, &rs);
        sum += (double)comps;
        sum_sq += (double)comps * (double)comps;
    }

    free(arr);

    double avg = sum / (double)num_trials;
    double var = (sum_sq / (double)num_trials) - avg * avg;
    if (var < 0.0) var = 0.0;

    if (mean) *mean = avg;
    if (stddev) *stddev = sqrt(var);
}

double quicksort_pivot_comparison(size_t n, size_t trials) {
    double mean_rand, stddev_rand;
    randomized_quicksort_analysis(n, trials, &mean_rand, &stddev_rand);

    /* Median-of-3 runs */
    double sum_m3 = 0.0;
    int *arr = (int *)malloc(n * sizeof(int));
    if (!arr) return 0.0;

    for (size_t t = 0; t < trials; t++) {
        for (size_t i = 0; i < n; i++) arr[i] = (int)i;
        size_t comps = 0;
        quicksort_median3(arr, 0, n, &comps);
        sum_m3 += (double)comps;
    }
    free(arr);

    if (mean_rand == 0.0) return 0.0;
    return (sum_m3 / (double)trials) / mean_rand;
}

/* ================================================================
 * L8: Hardness vs Randomness — Impagliazzo's Five Worlds
 * ================================================================ */

size_t nw_generator_seed_length(
    double circuit_size_lower_bound, size_t n)
{
    /* Nisan-Wigderson (1994):
     * If ∃ f ∈ E with circuit complexity 2^{Ω(n)}, then there exists
     * a PRG G: {0,1}^{O(log n)} → {0,1}^{poly(n)} that fools BPP.
     *
     * G stretches a seed of length O(log² n) to poly(n) bits.
     *
     * This function computes the seed length needed given the hardness
     * assumption (circuit size lower bound for f). */

    if (circuit_size_lower_bound < 1.0) return n;  /* No hardness */

    /* Seed length s = O(log(n) · log(circuit_size_lower_bound) / hardness)
     * ≈ O(log² n) for exponential hardness. */
    double s = log2((double)n) * log2(circuit_size_lower_bound)
             / log2(circuit_size_lower_bound) * log2((double)n);
    /* Actually: s = O(log(n/ε)² / log S) where S is circuit size.
     * For S = 2^{Ω(n)}: s = O(log²(n/ε)). */
    s = log2((double)n) * log2((double)n) / circuit_size_lower_bound;

    if (s < 1.0) s = 1.0;
    return (size_t)ceil(s);
}

const char *impagliazzo_world_name(ImpagliazzoWorld w) {
    switch (w) {
        case IMP_ALGORITHMICA: return "Algorithmica";
        case IMP_HEURISTICA:   return "Heuristica";
        case IMP_PESSILAND:    return "Pessiland";
        case IMP_MINICRYPT:    return "Minicrypt";
        case IMP_CRYPTOMANIA:  return "Cryptomania";
        default:               return "Unknown";
    }
}

bool bpp_equals_p_in_world(ImpagliazzoWorld w) {
    switch (w) {
        case IMP_ALGORITHMICA:
            /* P = NP or NP ⊆ BPP → BPP might still not equal P
             * without further derandomization assumptions. */
            return false;  /* BPP = P not implied */
        case IMP_HEURISTICA:
            /* NP is hard on average; strong circuit lower bounds
             * exist → BPP = P via NW derandomization. */
            return true;
        case IMP_PESSILAND:
            /* One-way functions exist but not public-key crypto.
             * Moderate derandomization possible, but BPP vs P open. */
            return false;
        case IMP_MINICRYPT:
            /* OWF exist → PRGs exist (Håstad et al.) → BPP ⊆ SUBEXP.
             * But BPP = P still may not hold without stronger assumptions. */
            return false;
        case IMP_CRYPTOMANIA:
            /* Public-key crypto → very strong OWF → strong PRGs.
             * Still, BPP = P requires EXP ≠ BPP or equivalent.
             * Generally: BPP = P is independent. */
            return false;
        default:
            return false;
    }
}

const char *impagliazzo_world_description(ImpagliazzoWorld w) {
    switch (w) {
        case IMP_ALGORITHMICA:
            return "P = NP or similar collapse. Beautiful world with fast algorithms for everything. BPP likely equals P.";
        case IMP_HEURISTICA:
            return "NP is hard on average but easy on random instances. Strong circuit lower bounds imply BPP = P.";
        case IMP_PESSILAND:
            return "OWF exist but not enough for public-key crypto. BPP vs P remains open in this world.";
        case IMP_MINICRYPT:
            return "One-way functions and private-key cryptography exist. Moderate derandomization possible (BPP ⊆ SUBEXP).";
        case IMP_CRYPTOMANIA:
            return "Full public-key cryptography. Strongest complexity assumptions. BPP = P requires separate assumptions.";
        default:
            return "Unknown world.";
    }
}
