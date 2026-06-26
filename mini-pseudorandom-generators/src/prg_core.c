/******************************************************************************
 * prg_core.c — Pseudorandom Generator Core
 *
 * Implements fundamental definitions, theorems, and constructions of
 * pseudorandom generators from complexity theory.
 *
 * Knowledge Coverage:
 *   L1: PRG, Distinguisher, NextBitPredictor definitions
 *   L2: Statistical distance, distinguish advantage, stretch, negligibility
 *   L3: Hybrid argument structure, polynomial stretch composition
 *   L4: Yao Theorem (1982), Hybrid Argument, IW97 derandomization
 *   L5: Statistical tests, stream encryption, cycle detection, NW generator
 *   L7: Stream cipher, BPP derandomization, cryptographic security
 *
 * References: Yao (FOCS 1982), Blum-Micali (1984), Goldreich (2001),
 *   Arora-Barak (2009), Impagliazzo-Wigderson (STOC 1997)
 * Courses: MIT 6.845, Stanford CS254, Princeton COS 522, Berkeley CS278
 ******************************************************************************/
#include "prg_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* =========================================================================
 * L1: Core Definitions - PRG formal model
 * PRG G:{0,1}^n -> {0,1}^{l(n)}, l(n)>n, PPT D: negligible advantage.
 * ========================================================================= */

double prg_statistical_distance(const double *p, const double *q, size_t n) {
    if (p == NULL || q == NULL || n == 0) return 0.0;
    double total = 0.0;
    for (size_t i = 0; i < n; i++) {
        double diff = p[i] - q[i];
        total += (diff < 0.0) ? -diff : diff;
    }
    return 0.5 * total;
}

double prg_distinguish_advantage(const bool *g_out, const bool *u_out,
                                  size_t n, size_t trials) {
    if (g_out == NULL || u_out == NULL || n == 0 || trials == 0) return 0.0;
    size_t g_ones = 0, u_ones = 0;
    for (size_t t = 0; t < trials; t++) {
        for (size_t i = 0; i < n; i++) {
            if (g_out[t * n + i]) g_ones++;
            if (u_out[t * n + i]) u_ones++;
        }
    }
    double denom = (double)(trials * n);
    if (denom == 0.0) return 0.0;
    double g_rate = (double)g_ones / denom;
    double u_rate = (double)u_ones / denom;
    double diff = g_rate - u_rate;
    return (diff < 0.0) ? -diff : diff;
}

/* =========================================================================
 * L4: Yao's Theorem (1982) - Next-bit unpredictability iff pseudorandomness
 * Theorem: {X_n} pseudorandom iff all next-bit tests pass with prob ~1/2.
 * ========================================================================= */

bool prg_yao_theorem_verify(const bool *prefix, size_t prefix_len,
                             bool predicted_bit, bool actual_bit,
                             size_t trials) {
    if (prefix == NULL || trials == 0) return false;
    size_t correct = 0;
    for (size_t t = 0; t < trials; t++) {
        if (predicted_bit == actual_bit) correct++;
    }
    double success_rate = (double)correct / (double)trials;
    double advantage = success_rate - 0.5;
    double threshold = 1.0 / sqrt((double)trials);
    return (advantage > threshold);
}

double prg_hybrid_argument_bound(size_t m, double nbp_adv) {
    if (m == 0) return 0.0;
    return (double)m * nbp_adv;
}

bool prg_has_stretch(const PRG *prg) {
    if (prg == NULL) return false;
    return prg->output_len > prg->seed_len;
}

double prg_stretch_ratio(const PRG *prg) {
    if (prg == NULL || prg->seed_len == 0) return 0.0;
    return (double)prg->output_len / (double)prg->seed_len;
}

/* =========================================================================
 * L5: Cycle Detection - Floyd's tortoise-and-hare algorithm
 * Finite-state PRGs eventually cycle. Floyd: O(1) space, O(mu+lambda) time.
 * ========================================================================= */

size_t prg_cycle_detect(const PRG *prg, const bool *seed, size_t *period) {
    if (prg == NULL || seed == NULL || prg->gen == NULL || period == NULL)
        return 0;
    size_t n = prg->seed_len;
    bool *slow = (bool *)malloc(n * sizeof(bool));
    bool *fast = (bool *)malloc(n * sizeof(bool));
    bool *temp = (bool *)malloc(n * sizeof(bool));
    if (slow == NULL || fast == NULL || temp == NULL) {
        free(slow); free(fast); free(temp);
        *period = 0; return 0;
    }
    memcpy(slow, seed, n * sizeof(bool));
    memcpy(fast, seed, n * sizeof(bool));
    size_t steps = 0;
    const size_t MAX_S = 1000000;
    bool found = false;
    while (steps < MAX_S) {
        prg->gen(slow, temp, n, n);
        memcpy(slow, temp, n * sizeof(bool));
        prg->gen(fast, temp, n, n);
        memcpy(fast, temp, n * sizeof(bool));
        prg->gen(fast, temp, n, n);
        memcpy(fast, temp, n * sizeof(bool));
        steps++;
        if (memcmp(slow, fast, n * sizeof(bool)) == 0) { found = true; break; }
    }
    if (found) {
        memcpy(slow, seed, n * sizeof(bool));
        size_t prefix_len = 0;
        while (prefix_len < MAX_S) {
            prg->gen(slow, temp, n, n);
            memcpy(slow, temp, n * sizeof(bool));
            prg->gen(fast, temp, n, n);
            memcpy(fast, temp, n * sizeof(bool));
            prefix_len++;
            if (memcmp(slow, fast, n * sizeof(bool)) == 0) break;
        }
        *period = 1;
        while (*period < MAX_S) {
            prg->gen(slow, temp, n, n);
            memcpy(slow, temp, n * sizeof(bool));
            if (memcmp(slow, fast, n * sizeof(bool)) == 0) break;
            (*period)++;
        }
    } else {
        *period = steps;
    }
    free(slow); free(fast); free(temp);
    return steps;
}
/* =========================================================================
 * L5: BPP Derandomization via PRG
 * Theorem: If PRG with seed O(log n) and poly output, then BPP subseteq SUBEXP.
 * IW97: If E has 2^{Omega(n)} circuits, then BPP = P.
 * ========================================================================= */

size_t prg_bpp_trials_for_error(double delta) {
    if (delta >= 1.0 || delta <= 0.0) return 1;
    double k = 18.0 * log(1.0 / delta);
    if (k < 1.0) k = 1.0;
    return (size_t)ceil(k);
}

void prg_stream_encrypt(const bool *key, size_t kl,
                         const bool *pt, size_t ml,
                         bool *ct, const PRG *prg) {
    if (key == NULL || pt == NULL || ct == NULL || prg == NULL) return;
    if (prg->gen == NULL) return;
    bool *ks = (bool *)malloc(ml * sizeof(bool));
    if (ks == NULL) return;
    prg->gen(key, ks, kl, ml);
    for (size_t i = 0; i < ml; i++) {
        ct[i] = pt[i] ^ ks[i];
    }
    free(ks);
}

bool prg_derandomize_bpp(bool (*dec)(const char*,size_t,const bool*,size_t),
                          const char *inp, size_t ilen, const PRG *prg) {
    if (dec == NULL || inp == NULL || prg == NULL) return false;
    if (prg->seed_len > 20) return false;
    size_t d = prg->seed_len, m = prg->output_len;
    size_t num_seeds = (size_t)1 << d;
    bool *seed = (bool *)calloc(d, sizeof(bool));
    bool *rb = (bool *)malloc(m * sizeof(bool));
    if (seed == NULL || rb == NULL) { free(seed); free(rb); return false; }
    size_t accepts = 0, rejects = 0;
    for (size_t s = 0; s < num_seeds; s++) {
        for (size_t i = 0; i < d; i++) seed[d-1-i] = (s>>i)&1;
        prg->gen(seed, rb, d, m);
        if (dec(inp, ilen, rb, m)) accepts++; else rejects++;
    }
    free(seed); free(rb);
    return accepts >= rejects;
}

/* =========================================================================
 * L2: Additional core concepts
 * ========================================================================= */

double prg_computational_distance(const bool *seq_a, const bool *seq_b,
                                   size_t n) {
    if (seq_a == NULL || seq_b == NULL || n == 0) return 0.0;
    size_t diff = 0;
    for (size_t i = 0; i < n; i++)
        if (seq_a[i] != seq_b[i]) diff++;
    return (double)diff / (double)n;
}

void prg_seed_expand(const PRG *prg, const bool *seed, bool *output) {
    if (prg == NULL || seed == NULL || output == NULL) return;
    if (prg->gen == NULL) return;
    prg->gen(seed, output, prg->seed_len, prg->output_len);
}

/* L3: Polynomial stretch from 1-bit stretch construction */

void prg_compose_stretch(const PRG *prg, const bool *seed, size_t k,
                          bool *output) {
    if (prg == NULL || seed == NULL || output == NULL || k == 0) return;
    if (prg->gen == NULL) return;
    size_t n = prg->seed_len;
    bool *state = (bool *)malloc(n * sizeof(bool));
    bool *next = (bool *)malloc((n+1) * sizeof(bool));
    if (state == NULL || next == NULL) { free(state); free(next); return; }
    memcpy(state, seed, n * sizeof(bool));
    for (size_t i = 0; i < k; i++) {
        prg->gen(state, next, n, n+1);
        output[i] = next[n];
        memcpy(state, next, n * sizeof(bool));
    }
    free(state); free(next);
}

/* L4: IW97 Framework */

double prg_iw97_bound(size_t m, double epsilon, size_t input_len) {
    if (epsilon <= 0.0 || m == 0) return 0.0;
    double seed_len = (double)m / epsilon;
    return pow(2.0, seed_len) * (double)input_len * (double)input_len;
}

/* =========================================================================
 * L7: Cryptographic PRG Security Applications
 * ========================================================================= */

double prg_security_level(size_t seed_len, double advantage) {
    if (advantage <= 0.0) return (double)seed_len;
    double t = 1.0 / advantage;
    if (t <= 1.0) return 1.0;
    return log2(t);
}

size_t prg_min_seed_length(size_t target_security, double advantage) {
    if (advantage <= 0.0) return target_security;
    double loss = -log2(advantage);
    if (loss < 0.0) loss = 0.0;
    return (size_t)ceil((double)target_security + loss);
}

bool prg_forward_secrecy_check(const bool *prev_output,
                                const bool *curr_state, size_t seed_len) {
    if (prev_output == NULL || curr_state == NULL || seed_len == 0) return false;
    size_t diff = 0;
    for (size_t i = 0; i < seed_len; i++)
        if (prev_output[i] != curr_state[i]) diff++;
    return diff > 0;
}

bool prg_is_negligible(double advantage, size_t n) {
    if (n == 0) return advantage < 1e-15;
    double threshold = pow(2.0, -(double)n);
    if (threshold < 1e-15) threshold = 1e-15;
    return advantage < threshold;
}

bool prg_compose(PRG *prgout, const PRG *g1, const PRG *g2) {
    if (prgout == NULL || g1 == NULL || g2 == NULL) return false;
    if (g1->output_len < g2->seed_len) return false;
    prgout->seed_len = g1->seed_len;
    prgout->output_len = g2->output_len;
    prgout->stretch = g2->output_len - g1->seed_len;
    prgout->adv = (g1->adv > g2->adv) ? g1->adv : g2->adv;
    return true;
}
size_t prg_bit_generation_rate(const PRG *prg) {
    if (prg == NULL) return 0;
    return prg->output_len;
}

size_t prg_entropy_loss(const PRG *prg) {
    if (prg == NULL) return 0;
    if (prg->output_len < prg->seed_len) return 0;
    return prg->output_len - prg->seed_len;
}

void prg_initialize(PRG *prg, size_t seed_len, size_t output_len,
                     void (*gen)(const bool*,bool*,size_t,size_t),
                     const char *name) {
    if (prg == NULL) return;
    prg->seed_len = seed_len;
    prg->output_len = output_len;
    prg->stretch = (output_len > seed_len) ? output_len - seed_len : 0;
    prg->gen = gen;
    prg->adv = 0.0;
    prg->name = name;
}

void prg_create_distinguisher(Distinguisher *d,
                               bool (*pred)(const bool*,size_t),
                               double adv) {
    if (d == NULL) return;
    d->pred = pred;
    d->queries = 0;
    d->adv = adv;
}

void prg_create_predictor(NextBitPredictor *nbp,
                           bool (*pred)(const bool*,size_t),
                           double prob, size_t pos) {
    if (nbp == NULL) return;
    nbp->pred_next = pred;
    nbp->prob = prob;
    nbp->pos = pos;
}

void prg_print_info(const PRG *prg) {
    if (prg == NULL) return;
    printf("PRG: %s\n", prg->name ? prg->name : "(unnamed)");
    printf("  Seed: %zu  Output: %zu  Stretch: %zu  Ratio: %.2f\n",
           prg->seed_len, prg->output_len, prg->stretch,
           prg_stretch_ratio(prg));
    printf("  Advantage: %.2e\n", prg->adv);
}

/* =========================================================================
 * L5: Nisan-Wigderson Generator Framework
 * NW generator: {0,1}^d -> {0,1}^m based on hard function f in E
 * and combinatorial design S_1,...,S_m subseteq [d] with small intersections.
 * NW(x)_i = f(x|_{S_i}). Key tool behind IW97.
 * ========================================================================= */

size_t prg_nw_design_size(size_t m, size_t k, size_t r) {
    if (r == 0 || k == 0) return m;
    double d_min = pow((double)m, 1.0/(double)(r+1));
    size_t d = (size_t)ceil(d_min);
    return (d < k) ? k : d;
}

size_t prg_nw_output_length(size_t d, size_t r) {
    if (d == 0) return 0;
    size_t m = 1;
    for (size_t i = 0; i <= r && m <= SIZE_MAX/d; i++) m *= d;
    return m;
}

size_t prg_nw_hardness_to_stretch(size_t hardness, size_t seed_length) {
    if (hardness <= seed_length) return seed_length + 1;
    size_t gap = hardness - seed_length;
    size_t extra = (size_t)1 << (gap < 32 ? gap : 31);
    return seed_length + extra;
}

double prg_security_parameter(size_t seed_len, size_t time_bound) {
    double s = (double)seed_len, t = (double)time_bound;
    return (s < t) ? s : t;
}

double prg_advantage_decay(double initial_adv, size_t trials) {
    if (trials == 0) return initial_adv;
    double decay = exp(-2.0 * (double)trials * initial_adv * initial_adv);
    double result = initial_adv * decay;
    return (result < 0.0) ? 0.0 : result;
}

double prg_unpredictability_to_indistinguishability(double nbp_advantage,
                                                      size_t num_bits) {
    return prg_hybrid_argument_bound(num_bits, nbp_advantage);
}

double prg_indistinguishability_to_unpredictability(double distinguish_adv,
                                                       size_t num_bits) {
    if (num_bits == 0) return distinguish_adv;
    return distinguish_adv / (double)num_bits;
}

void prg_quantitative_yao(double nbp_adv, double dist_adv, size_t num_bits,
                           double *out_nbp_bound, double *out_dist_bound) {
    if (out_nbp_bound)
        *out_nbp_bound = prg_unpredictability_to_indistinguishability(nbp_adv, num_bits);
    if (out_dist_bound)
        *out_dist_bound = prg_indistinguishability_to_unpredictability(dist_adv, num_bits);
}

/* =========================================================================
 * Extended PRG Analysis Functions
 * ========================================================================= */

/**
 * prg_entropy_estimate - Estimate min-entropy from bit sequence.
 *
 * Uses the empirical distribution of t-bit blocks.
 * H_inf = -log2(max_{block} Pr[block]).
 */
double prg_entropy_estimate(const bool *bits, size_t n, size_t block_size) {
    if (bits == NULL || n == 0 || block_size == 0 || block_size > 8) return 0.0;
    size_t num_blocks = n / block_size;
    if (num_blocks == 0) return 0.0;
    size_t num_patterns = (size_t)1 << block_size;
    size_t *counts = (size_t *)calloc(num_patterns, sizeof(size_t));
    if (counts == NULL) return 0.0;
    for (size_t i = 0; i < num_blocks; i++) {
        size_t idx = 0;
        for (size_t j = 0; j < block_size; j++)
            idx = (idx << 1) | (bits[i * block_size + j] ? 1 : 0);
        counts[idx]++;
    }
    size_t max_count = 0;
    for (size_t i = 0; i < num_patterns; i++)
        if (counts[i] > max_count) max_count = counts[i];
    free(counts);
    if (max_count == 0) return 0.0;
    double max_prob = (double)max_count / (double)num_blocks;
    return -log2(max_prob);
}

/**
 * prg_compression_ratio - Estimate compressibility of PRG output.
 *
 * Compression ratio = 1 - entropy_rate. Lower is better (harder to compress).
 * Encodes the sequence using simple run-length encoding and compares.
 */
double prg_compression_ratio(const bool *bits, size_t n) {
    if (bits == NULL || n == 0) return 0.0;
    /* Simple RLE: count runs */
    size_t runs = 1;
    for (size_t i = 1; i < n; i++)
        if (bits[i] != bits[i-1]) runs++;
    /* RLE compressed size: runs * (1 + log2(max_run_length)) bits */
    double rle_bits = (double)runs * (1.0 + log2((double)n));
    return rle_bits / (double)n;
}

/**
 * prg_correlation_matrix - Compute pairwise correlation between bit positions.
 *
 * For positions (i,j), computes the correlation of b_i and b_j across trials.
 * Returns the maximum absolute correlation found.
 */
double prg_correlation_matrix(const bool *bits, size_t n, size_t num_positions) {
    if (bits == NULL || n == 0 || num_positions < 2) return 0.0;
    if (num_positions > 64) num_positions = 64;
    /* Treat the bits as consecutive blocks of length num_positions */
    size_t num_samples = n / num_positions;
    if (num_samples < 2) return 0.0;
    double max_corr = 0.0;
    for (size_t i = 0; i < num_positions; i++) {
        for (size_t j = i + 1; j < num_positions; j++) {
            double sum_ij = 0.0, sum_i = 0.0, sum_j = 0.0;
            double sum_ii = 0.0, sum_jj = 0.0;
            for (size_t s = 0; s < num_samples; s++) {
                double bi = bits[s * num_positions + i] ? 1.0 : 0.0;
                double bj = bits[s * num_positions + j] ? 1.0 : 0.0;
                sum_ij += bi * bj;
                sum_i += bi; sum_j += bj;
                sum_ii += bi * bi; sum_jj += bj * bj;
            }
            double num = (double)num_samples * sum_ij - sum_i * sum_j;
            double den1 = (double)num_samples * sum_ii - sum_i * sum_i;
            double den2 = (double)num_samples * sum_jj - sum_j * sum_j;
            double den = sqrt(den1 * den2);
            double corr = (den > 0) ? fabs(num / den) : 0.0;
            if (corr > max_corr) max_corr = corr;
        }
    }
    return max_corr;
}

/**
 * prg_mutual_information - Approximate mutual information between
 * first half and second half of the PRG output.
 *
 * I(X;Y) = H(X) + H(Y) - H(X,Y)
 * If X and Y are independent, I(X;Y) = 0.
 */
double prg_mutual_information(const bool *bits, size_t n) {
    if (bits == NULL || n < 2) return 0.0;
    size_t half = n / 2;
    size_t ones1 = 0, ones2 = 0, joint_11 = 0;
    for (size_t i = 0; i < half; i++) {
        bool b1 = bits[i];
        bool b2 = bits[half + i];
        if (b1) ones1++;
        if (b2) ones2++;
        if (b1 && b2) joint_11++;
    }
    double p1 = (double)ones1 / (double)half;
    double p2 = (double)ones2 / (double)half;
    double p11 = (double)joint_11 / (double)half;
    double p10 = p1 - p11, p01 = p2 - p11, p00 = 1.0 - p1 - p2 + p11;
    double hx = 0.0, hy = 0.0, hxy = 0.0;
    if (p1 > 0 && p1 < 1) hx = -p1 * log2(p1) - (1-p1) * log2(1-p1);
    if (p2 > 0 && p2 < 1) hy = -p2 * log2(p2) - (1-p2) * log2(1-p2);
    if (p00 > 0) hxy -= p00 * log2(p00);
    if (p01 > 0) hxy -= p01 * log2(p01);
    if (p10 > 0) hxy -= p10 * log2(p10);
    if (p11 > 0) hxy -= p11 * log2(p11);
    double mi = hx + hy - hxy;
    return (mi < 0) ? 0.0 : mi;
}

/**
 * prg_chi_squared_uniformity - Chi-squared test for output uniformity.
 *
 * Divides output into k buckets and tests if each bucket is equally likely.
 * chi^2 = sum (O_i - E)^2 / E where E = n/k.
 */
double prg_chi_squared_uniformity(const bool *bits, size_t n, size_t buckets) {
    if (bits == NULL || n == 0 || buckets == 0) return 0.0;
    size_t bits_per_bucket = 8;
    size_t vals_per_bucket = 0;
    if (buckets <= 256) {
        bits_per_bucket = 8;
        vals_per_bucket = 256 / buckets;
        if (vals_per_bucket == 0) vals_per_bucket = 1;
    }
    size_t *counts = (size_t *)calloc(buckets, sizeof(size_t));
    if (counts == NULL) return 0.0;
    for (size_t i = 0; i + bits_per_bucket <= n; i += bits_per_bucket) {
        size_t val = 0;
        for (size_t j = 0; j < bits_per_bucket; j++)
            val = (val << 1) | (bits[i+j] ? 1 : 0);
        size_t bucket = val / vals_per_bucket;
        if (bucket < buckets) counts[bucket]++;
    }
    size_t total = n / bits_per_bucket;
    double expected = (double)total / (double)buckets;
    double chi_sq = 0.0;
    for (size_t i = 0; i < buckets; i++) {
        double diff = (double)counts[i] - expected;
        if (expected > 0) chi_sq += diff * diff / expected;
    }
    free(counts);
    return chi_sq;
}

/**
 * prg_poker_test - Poker test for PRG output randomness.
 *
 * Classic poker test: counts frequencies of m-bit blocks and
 * compares to the distribution expected for random data.
 * X = (2^m / k) * sum(f_i^2) - k where k = floor(n/m).
 */
double prg_poker_test(const bool *bits, size_t n, size_t m) {
    if (bits == NULL || n == 0 || m == 0 || m > 8) return 0.0;
    size_t k = n / m;
    if (k == 0) return 0.0;
    size_t num_patterns = (size_t)1 << m;
    size_t *counts = (size_t *)calloc(num_patterns, sizeof(size_t));
    if (counts == NULL) return 0.0;
    for (size_t i = 0; i < k; i++) {
        size_t idx = 0;
        for (size_t j = 0; j < m; j++)
            idx = (idx << 1) | (bits[i * m + j] ? 1 : 0);
        counts[idx]++;
    }
    double sum_sq = 0.0;
    for (size_t i = 0; i < num_patterns; i++)
        sum_sq += (double)counts[i] * (double)counts[i];
    free(counts);
    double x = ((double)num_patterns / (double)k) * sum_sq - (double)k;
    return x;
}
