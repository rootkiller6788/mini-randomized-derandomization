/**
 * @file nw_prg.c
 * @brief Nisan-Wigderson PRG: Construction, evaluation, derandomization, analysis.
 *
 * Implements:
 *   L2: PRG construction strategies (standard, XOR-amplified, IW97, optimized)
 *   L4: Stretch vs hardness tradeoff computation
 *   L5: PRG evaluation (full, lazy, word-based)
 *   L5: Distinguisher advantage estimation
 *   L5: Hybrid argument for the NW proof
 *   L7: BPP derandomization via NW PRG
 *   L7: Cryptographic applications (key derivation)
 *   L5: Statistical testing suite for PRG output
 *
 * Reference:
 *   Nisan & Wigderson, "Hardness vs Randomness", JCSS 49(2), 1994
 *   Impagliazzo & Wigderson, "P = BPP if E requires exponential circuits", STOC 1997
 *   Arora & Barak, Ch 20
 *   Goldreich, "Computational Complexity: A Conceptual Perspective", Ch 8
 */

#include "nw_prg.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

/* =============================================
 * Default hard function: Inner Product mod 2
 *
 * IP(x_1,...,x_n) = XOR_{i=1}^{n/2} (x_i AND x_{i+n/2})
 *
 * This is the Goldreich-Levin hard-core predicate: if f is a
 * one-way function, IP(x, r) is unpredictable given f(x), r.
 *
 * For the NW construction, this serves as a candidate hard
 * function under standard cryptographic assumptions.
 * ============================================= */

bool nw_default_hard_function(const bool *input, size_t n) {
    size_t half = n / 2;
    bool result = false;
    for (size_t i = 0; i < half; i++) {
        result ^= (input[i] && input[i + half]);
    }
    return result;
}

/* =============================================
 * L4: Optimal PRG Construction
 *
 * Given desired output length m and hardness guarantee S,
 * construct the optimal NW PRG using:
 *   1. Parameter optimization (k, l, intersect_bound)
 *   2. Design construction (Reed-Solomon when possible, random otherwise)
 *   3. Optional hardness amplification (Yao XOR)
 * ============================================= */

NWPRG *nw_prg_construct_optimal(size_t output_len, double hardness, NWPRGMode mode) {
    if (output_len == 0 || hardness <= 0) return NULL;

    /* Step 1: Determine design parameters */
    size_t m = output_len;
    if (m < 2) m = 2;

    double log_S = log2(hardness);
    if (log_S < 1.0) log_S = 1.0;
    size_t l = (size_t)ceil(log_S);
    if (l < 2) l = 2;

    /* Intersection bound = ceil(log m) */
    size_t intersect_bound = nw_ceil_log2(m);
    if (intersect_bound < 1) intersect_bound = 1;

    /* Seed length k = l^2 (for Reed-Solomon designs) */
    size_t k = l * l;
    size_t log_m = nw_ceil_log2(m);
    if (k < log_m + l) k = log_m + l;

    /* Step 2: Construct design */
    NWDesign *design = NULL;

    switch (mode) {
        case NW_PRG_STANDARD:
            /* Try Reed-Solomon first */
            design = nw_reed_solomon_design_param(k, m, l, intersect_bound);
            if (!design) {
                /* Fall back to random design */
                design = nw_random_design(k, m, l);
            }
            break;

        case NW_PRG_YAO_XOR:
            /* First amplify hardness, then construct */
            {
                /* XOR-k to boost hardness: choose k s.t. epsilon^k <= 1/m */
                double epsilon = 0.1; /* base epsilon */
                size_t xor_k = 1;
                while (pow(epsilon, (double)xor_k) > 1.0 / (double)m && xor_k < 10)
                    xor_k++;
                /* Amortize: adjust l to be l * xor_k effectively */
                design = nw_reed_solomon_design_param(k, m, l * xor_k, intersect_bound * xor_k);
                if (!design)
                    design = nw_random_design(k, m, l);
            }
            break;

        case NW_PRG_IW97:
            /* IW97: Full derandomization pipeline
             * Requires exponential hardness: S = 2^{Omega(n)} */
            {
                size_t k_iw = (size_t)ceil(log((double)m) * log((double)m));
                if (k_iw < 2) k_iw = 2;
                design = nw_random_design(k_iw, m, k_iw / 2);
                if (!design) design = nw_random_design(k, m, l);
                k = k_iw;
            }
            break;

        case NW_PRG_OPTIMIZED:
            /* Optimized: balance all parameters for maximum stretch */
            {
                /* Grid search over feasible (k,l) pairs */
                size_t best_k = k;
                double best_stretch = 0;
                for (size_t try_l = 2; try_l <= l + 10 && try_l <= 256; try_l++) {
                    size_t try_k = try_l * try_l;
                    size_t try_m = nw_output_length(try_k, hardness);
                    double stretch = (double)try_m / (double)try_k;
                    if (stretch > best_stretch) {
                        best_stretch = stretch;
                        best_k = try_k;
                    }
                }
                k = best_k;
                design = nw_reed_solomon_design_param(k, m, l, intersect_bound);
                if (!design) design = nw_random_design(k, m, l);
            }
            break;

        default:
            design = nw_random_design(k, m, l);
    }

    if (!design) return NULL;

    /* Step 3: Create hard function (Goldreich-Levin inner product predicate) */
    /* For a real implementation, this would be e.g.:
     *   - Goldreich-Levin hard-core bit
     *   - Subset sum hard function
     *   - Discrete log hard predicate
     */
    /* Default hard function: inner product mod 2 (provably hard under assumptions) */
    HardFunction *hf = nw_hard_function_create(
        design->l, nw_default_hard_function,
        hardness, 0.1, "inner-product-hard-function"
    );

    if (!hf) { nw_design_free(design); return NULL; }

    /* Step 4: Assemble PRG */
    NWPRG *prg = nw_prg_construct(design, hf);
    if (!prg) {
        nw_hard_function_free(hf);
        nw_design_free(design);
        return NULL;
    }

    return prg;
}

/* =============================================
 * L4: Stretch / Seed Trade-off Computation
 *
 * The stretch m/k is the key efficiency metric:
 *   - Standard NW: m/k = 2^{Theta(sqrt(log S))} = quasipolynomial
 *   - Exponential hardness: m/k = 2^{Omega(k)} (exponential stretch)
 * ============================================= */

size_t nw_prg_max_stretch(size_t seed_len, double hardness) {
    if (seed_len == 0 || hardness <= 0) return 0;

    /* log hardness controls the design parameters.
     * With seed_len = k, we have l ≈ sqrt(k), and
     * the design supports m ≈ l^{l} ≈ 2^{sqrt(k) log k} sets. */

    double log_S = log2(hardness);
    if (log_S < 1.0) log_S = 1.0;

    /* Check if hardness is exponential: S = 2^{Omega(n)} */
    bool is_exp = (log_S > 2.0 * (double)seed_len);

    if (is_exp) {
        /* Exponential stretch: m = 2^{Omega(k)} */
        size_t stretch = (size_t)1 << (seed_len / 4);
        if (stretch > (size_t)1 << 30) stretch = (size_t)1 << 30;
        return stretch;
    } else {
        /* Sub-exponential: m = 2^{k^{1/3}} */
        double exp_part = pow((double)seed_len, 1.0 / 3.0);
        size_t stretch = (size_t)pow(2.0, exp_part);
        if (stretch < 2) stretch = 2;
        return stretch;
    }
}

size_t nw_prg_min_seed(size_t output_len, double hardness) {
    return nw_seed_length(output_len, hardness);
}

/* =============================================
 * L5: PRG Evaluation (Full, Range, Word-based)
 * ============================================= */

void nw_prg_eval_all(const NWPRG *prg, const bool *seed, bool *out) {
    nw_prg_evaluate(prg, seed, out);
}

void nw_prg_eval_range(const NWPRG *prg, const bool *seed,
                        size_t start_idx, size_t count, bool *out)
{
    if (!prg || !seed || !out) return;
    if (start_idx >= prg->output_len) return;
    if (start_idx + count > prg->output_len)
        count = prg->output_len - start_idx;

    bool *restricted = (bool *)malloc(prg->design->l * sizeof(bool));
    if (!restricted) return;

    for (size_t i = 0; i < count; i++) {
        size_t global_idx = start_idx + i;
        for (size_t j = 0; j < prg->design->l; j++) {
            size_t idx = prg->design->sets[global_idx][j];
            restricted[j] = (idx < prg->seed_len) ? seed[idx] : false;
        }
        out[i] = prg->hard_fn->evaluate(restricted, prg->design->l);
    }

    free(restricted);
}

void nw_prg_eval_word(const NWPRG *prg, uint64_t seed_word, bool *out) {
    if (!prg || !out) return;
    bool *seed = (bool *)calloc(prg->seed_len, sizeof(bool));
    if (!seed) return;
    nw_int_to_seed(prg->seed_len, seed_word, seed);
    nw_prg_evaluate(prg, seed, out);
    free(seed);
}

uint64_t nw_prg_eval_word_to_word(const NWPRG *prg, uint64_t seed_word) {
    if (!prg || prg->output_len > 64) return 0;
    bool *out = (bool *)calloc(prg->output_len, sizeof(bool));
    if (!out) return 0;
    nw_prg_eval_word(prg, seed_word, out);
    uint64_t result = 0;
    for (size_t i = 0; i < prg->output_len && i < 64; i++)
        if (out[i]) result |= (1ULL << i);
    free(out);
    return result;
}

/* =============================================
 * L5: Distinguisher Advantage Estimation
 *
 * Given a distinguisher circuit D: {0,1}^m -> {0,1},
 * estimate its advantage:
 *   ε = |Pr_{x<-U_k}[D(G(x))=1] - Pr_{y<-U_m}[D(y)=1]|
 *
 * Uses Monte Carlo sampling for each distribution.
 * ============================================= */

double nw_prg_distinguisher_advantage(
    const NWPRG *prg,
    const Circuit *distinguish_circuit,
    size_t num_samples)
{
    if (!prg || !distinguish_circuit || num_samples == 0) return 0.0;

    size_t count_prg = 0;
    size_t count_random = 0;

    bool *seed = (bool *)malloc(prg->seed_len * sizeof(bool));
    bool *prg_out = (bool *)malloc(prg->output_len * sizeof(bool));
    bool *random_out = (bool *)malloc(prg->output_len * sizeof(bool));
    bool output_val;

    if (!seed || !prg_out || !random_out) {
        free(seed); free(prg_out); free(random_out);
        return 0.0;
    }

    for (size_t s = 0; s < num_samples; s++) {
        /* PRG distribution */
        nw_random_bits(seed, prg->seed_len);
        nw_prg_evaluate(prg, seed, prg_out);
        nw_circuit_evaluate(distinguish_circuit, prg_out, &output_val);
        if (output_val) count_prg++;

        /* Random distribution */
        nw_random_bits(random_out, prg->output_len);
        nw_circuit_evaluate(distinguish_circuit, random_out, &output_val);
        if (output_val) count_random++;
    }

    free(seed); free(prg_out); free(random_out);

    double prob_prg = (double)count_prg / (double)num_samples;
    double prob_random = (double)count_random / (double)num_samples;

    return fabs(prob_prg - prob_random);
}

/* =============================================
 * L5: Distinguisher to Predictor (Yao)
 *
 * Yao (1982): If D distinguishes G(U_k) from U_m with
 * advantage ε, then for some i there is a next-bit predictor
 * with advantage ε/m.
 *
 * Construction:
 *   On input (y_1,...,y_i), P_i guesses:
 *     If D(y_1,...,y_i, r_{i+1},...,r_m) = 1, output r_{i+1}
 *     else output 1 - r_{i+1}
 *   where r's are random bits.
 *
 * By the hybrid argument, P_i's advantage is related to
 * the gap between hybrid H_i and H_{i-1}.
 * ============================================= */

NextBitPredictor *nw_distinguisher_to_predictors(
    const NWPRG *prg,
    const Circuit *distinguisher,
    double advantage)
{
    if (!prg || !distinguisher) return NULL;

    /* Find the hybrid position with largest gap */
    int best_idx = nw_hybrid_position(prg, distinguisher, 500);
    if (best_idx < 0) return NULL;

    NextBitPredictor *predictor = (NextBitPredictor *)malloc(sizeof(NextBitPredictor));
    if (!predictor) return NULL;

    predictor->n = prg->output_len;
    predictor->position = (size_t)best_idx;
    predictor->predictor = NULL; /* The predictor circuit would be built here */
    predictor->advantage = advantage / (double)prg->output_len;

    return predictor;
}

/**
 * L5: Hybrid Argument Position Search
 *
 * Hybrid H_i: first i bits from PRG, rest random.
 * Measure D on each hybrid, find maximum gap between adjacent hybrids.
 */
int nw_hybrid_position(const NWPRG *prg, const Circuit *distinguisher, size_t num_samples) {
    if (!prg || !distinguisher) return -1;

    size_t m = prg->output_len;
    if (m == 0) return -1;

    double *hybrid_probs = (double *)calloc(m + 1, sizeof(double));
    bool *hybrid_out = (bool *)malloc(m * sizeof(bool));
    bool *seed = (bool *)malloc(prg->seed_len * sizeof(bool));

    if (!hybrid_probs || !hybrid_out || !seed) {
        free(hybrid_probs); free(hybrid_out); free(seed);
        return -1;
    }

    for (size_t i = 0; i <= m; i++) {
        size_t count = 0;
        for (size_t s = 0; s < num_samples; s++) {
            /* First i bits from PRG */
            nw_random_bits(seed, prg->seed_len);
            nw_prg_eval_range(prg, seed, 0, i, hybrid_out);
            /* Remaining m-i bits random */
            nw_random_bits(hybrid_out + i, m - i);

            bool result;
            nw_circuit_evaluate(distinguisher, hybrid_out, &result);
            if (result) count++;
        }
        hybrid_probs[i] = (double)count / (double)num_samples;
    }

    /* Find largest gap */
    int best_idx = -1;
    double max_gap = 0.0;
    for (size_t i = 1; i <= m; i++) {
        double gap = fabs(hybrid_probs[i] - hybrid_probs[i-1]);
        if (gap > max_gap) {
            max_gap = gap;
            best_idx = (int)(i - 1);
        }
    }

    free(hybrid_probs); free(hybrid_out); free(seed);
    return best_idx;
}

/* =============================================
 * L7: BPP Derandomization via NW PRG
 *
 * Given a BPP algorithm A(x, r) using m random bits,
 * enumerate all 2^k seeds and take majority vote.
 *
 * Error analysis:
 *   If A has error <= 1/3 for each random string,
 *   and G fools circuits of appropriate size,
 *   then majority over PRG outputs is correct with
 *   probability >= 1 - exp(-2^k).
 * ============================================= */

bool nw_bpp_derandomize(
    const NWPRG *prg,
    bool (*bp_algorithm)(const bool *input, const bool *random, size_t n, size_t m),
    const bool *input, size_t n, size_t m,
    bool *result)
{
    if (!prg || !bp_algorithm || !input || !result) return false;

    size_t k = prg->seed_len;
    size_t num_seeds = (size_t)1 << k;

    /* Limit seeds for practical demonstration */
    if (num_seeds > (size_t)1 << 20) num_seeds = (size_t)1 << 16;

    bool *prg_output = (bool *)malloc(m * sizeof(bool));
    bool *seed = (bool *)malloc(k * sizeof(bool));
    if (!prg_output || !seed) {
        free(prg_output); free(seed);
        return false;
    }

    size_t votes_yes = 0;
    size_t votes_no = 0;

    for (size_t s = 0; s < num_seeds; s++) {
        nw_int_to_seed(k, s, seed);
        nw_prg_evaluate(prg, seed, prg_output);

        bool vote = bp_algorithm(input, prg_output, n, m);
        if (vote) votes_yes++;
        else votes_no++;
    }

    free(prg_output); free(seed);

    *result = (votes_yes > votes_no);
    return true;
}

NWPRG *nw_prg_for_BPP(size_t n, size_t randomness, double hardness) {
    (void)n;
    if (randomness == 0) return NULL;

    /* For a BPP algorithm using R random bits on inputs of size n:
     * We need a PRG stretching O(log n) bits to R bits.
     * Seed length k = O(log n) for P-level simulation,
     * or k = n^{o(1)} for SUBEXP-level simulation. */

    size_t k = nw_seed_length(randomness, hardness);
    if (k < 2) k = 2;

    NWPRG *prg = nw_prg_construct_optimal(randomness, hardness, NW_PRG_STANDARD);
    return prg;
}

/* =============================================
 * L5: Statistical Testing Suite
 *
 * For demonstration: apply basic statistical tests
 * to PRG output to verify pseudorandomness empirically.
 * (Note: passing statistical tests != cryptographic security)
 * ============================================= */

/* Frequency (monobit) test: count proportion of 1s */
static double _monobit_test(const bool *stream, size_t len) {
    size_t ones = 0;
    for (size_t i = 0; i < len; i++)
        if (stream[i]) ones++;
    double prop = (double)ones / (double)len;
    /* Test statistic: |prop - 0.5| */
    return fabs(prop - 0.5);
}

/* Runs test: count number of runs */
static double _runs_test(const bool *stream, size_t len) {
    if (len < 2) return 0.0;
    size_t runs = 1;
    for (size_t i = 1; i < len; i++)
        if (stream[i] != stream[i-1]) runs++;
    /* Expected runs = (2*len*ones*zeros)/len + 1 */
    size_t ones = 0;
    for (size_t i = 0; i < len; i++)
        if (stream[i]) ones++;
    size_t zeros = len - ones;
    double expected = 2.0 * (double)ones * (double)zeros / (double)len + 1.0;
    if (expected == 0) return 0.0;
    return fabs((double)runs - expected) / sqrt(expected);
}

/* Longest run test */
static double _longest_run_test(const bool *stream, size_t len) {
    size_t longest = 0, current = 0;
    bool prev = false;
    for (size_t i = 0; i < len; i++) {
        if (stream[i] == prev) {
            current++;
        } else {
            if (current > longest) longest = current;
            current = 1;
            prev = stream[i];
        }
    }
    if (current > longest) longest = current;
    /* For truly random bits, longest run ~ log2(n) - log2(log2(n)) */
    double expected_log = log2((double)len);
    return fabs((double)longest - expected_log) / expected_log;
}

StatisticalTestResult *nw_prg_statistical_tests(
    const NWPRG *prg, size_t num_streams, size_t stream_len)
{
    if (!prg || num_streams == 0 || stream_len == 0) return NULL;

    StatisticalTestResult *r = (StatisticalTestResult *)malloc(sizeof(StatisticalTestResult));
    if (!r) return NULL;

    r->num_samples = num_streams;
    r->num_tests = 3; /* monobit, runs, longest-run */
    r->p_values = (double *)calloc(r->num_tests, sizeof(double));
    r->passes_all = true;
    r->min_p_value = 1.0;

    if (!r->p_values) { free(r); return NULL; }

    bool *stream = (bool *)malloc(stream_len * sizeof(bool));
    bool *seed = (bool *)malloc(prg->seed_len * sizeof(bool));

    if (!stream || !seed) {
        free(stream); free(seed); free(r->p_values); free(r);
        return NULL;
    }

    double *test_stats[3];
    test_stats[0] = (double *)calloc(num_streams, sizeof(double));
    test_stats[1] = (double *)calloc(num_streams, sizeof(double));
    test_stats[2] = (double *)calloc(num_streams, sizeof(double));

    for (size_t s = 0; s < num_streams; s++) {
        /* Generate seed */
        nw_int_to_seed(prg->seed_len, s, seed);

        /* Fill stream by evaluating PRG repeatedly */
        size_t pos = 0;
        while (pos < stream_len) {
            size_t chunk = prg->output_len;
            if (pos + chunk > stream_len) chunk = stream_len - pos;
            nw_prg_eval_range(prg, seed, 0, chunk, stream + pos);
            pos += chunk;
            /* Increment seed for next chunk */
            uint64_t seed_val = nw_seed_to_int(prg->seed_len, seed);
            seed_val++;
            nw_int_to_seed(prg->seed_len, seed_val, seed);
        }

        test_stats[0][s] = _monobit_test(stream, stream_len);
        test_stats[1][s] = _runs_test(stream, stream_len);
        test_stats[2][s] = _longest_run_test(stream, stream_len);
    }

    /* Aggregate: average test statistic */
    for (size_t t = 0; t < 3; t++) {
        double avg = 0.0;
        for (size_t s = 0; s < num_streams; s++)
            avg += test_stats[t][s];
        avg /= (double)num_streams;
        r->p_values[t] = avg;
        if (avg < r->min_p_value) r->min_p_value = avg;
    }

    /* Pass if all statistics < threshold (simplified) */
    r->passes_all = (r->p_values[0] < 0.1 && r->p_values[1] < 0.5 && r->p_values[2] < 0.3);

    free(stream); free(seed);
    free(test_stats[0]); free(test_stats[1]); free(test_stats[2]);

    return r;
}

void nw_statistical_test_result_free(StatisticalTestResult *r) {
    if (!r) return;
    free(r->p_values);
    free(r);
}

/* =============================================
 * L7: Cryptographic Applications
 * ============================================= */

/**
 * Key Derivation Function (KDF):
 * Given a short random seed (key), expand to a keystream
 * of arbitrary length using the NW PRG in counter mode.
 *
 * This demonstrates the PRG's use as a stream cipher.
 * Security: if the hard function is truly hard, the
 * output is computationally indistinguishable from random.
 */
void nw_prg_key_derivation(
    const NWPRG *prg, const bool *seed,
    bool *keystream, size_t stream_len)
{
    if (!prg || !seed || !keystream) return;

    size_t k = prg->seed_len;
    bool *current_seed = (bool *)malloc(k * sizeof(bool));
    bool *chunk = (bool *)malloc(prg->output_len * sizeof(bool));

    if (!current_seed || !chunk) {
        free(current_seed); free(chunk);
        return;
    }

    memcpy(current_seed, seed, k * sizeof(bool));

    size_t pos = 0;
    while (pos < stream_len) {
        nw_prg_evaluate(prg, current_seed, chunk);
        size_t copy_len = prg->output_len;
        if (pos + copy_len > stream_len) copy_len = stream_len - pos;
        memcpy(keystream + pos, chunk, copy_len * sizeof(bool));
        pos += copy_len;

        /* Increment seed (counter mode) */
        uint64_t val = nw_seed_to_int(k, current_seed);
        val++;
        nw_int_to_seed(k, val, current_seed);
    }

    free(current_seed); free(chunk);
}

/**
 * Simulate a randomized algorithm using NW PRG:
 * Enumerate seeds, run algorithm with PRG output as randomness,
 * take majority vote.
 */
bool nw_prg_simulate_randomized(
    const NWPRG *prg,
    bool (*algorithm)(const void *input, const bool *random, size_t in_sz, size_t r_sz),
    const void *input, size_t input_size,
    size_t rand_bits, double error_bound)
{
    (void)error_bound;
    if (!prg || !algorithm || !input || rand_bits == 0) return false;

    size_t k = prg->seed_len;
    size_t num_seeds = (size_t)1 << k;
    if (num_seeds > (size_t)1 << 16) num_seeds = (size_t)1 << 16;

    bool *random_bits = (bool *)malloc(rand_bits * sizeof(bool));
    if (!random_bits) return false;

    size_t yes = 0, no = 0;

    for (size_t s = 0; s < num_seeds; s++) {
        /* Generate PRG output as randomness */
        size_t pos = 0;
        uint64_t seed_val = s;
        bool seed_bits[64];
        nw_int_to_seed(k, seed_val, seed_bits);

        while (pos < rand_bits) {
            size_t chunk = prg->output_len;
            if (pos + chunk > rand_bits) chunk = rand_bits - pos;
            /* Evaluate PRG for this seed */
            bool out[1024];
            nw_prg_eval_range(prg, seed_bits, 0, chunk, out);
            for (size_t c = 0; c < chunk; c++)
                random_bits[pos + c] = out[c];
            pos += chunk;
            seed_val++;
            nw_int_to_seed(k, seed_val, seed_bits);
        }

        bool vote = algorithm(input, random_bits, input_size, rand_bits);
        if (vote) yes++; else no++;
    }

    free(random_bits);
    return (yes > no);
}

/* =============================================
 * Utility: PRG Management
 * ============================================= */

size_t nw_prg_security_level(const NWPRG *prg) {
    if (!prg || !prg->hard_fn) return 0;
    return (size_t)log2(prg->hard_fn->hardness);
}

void nw_prg_print(const NWPRG *prg, FILE *fp) {
    if (!prg) { fprintf(fp, "NWPRG(NULL)\n"); return; }
    fprintf(fp, "NW-PRG {\n");
    fprintf(fp, "  seed_len=%zu, output_len=%zu, stretch=%.2f\n",
            prg->seed_len, prg->output_len, prg->stretch);
    fprintf(fp, "  design: k=%zu, m=%zu, l=%zu, intersect_bound=%zu\n",
            prg->design->k, prg->design->m, prg->design->l,
            prg->design->intersect_bound);
    if (prg->hard_fn) {
        fprintf(fp, "  hard_fn: n=%zu, hardness=%.1f, epsilon=%.4f\n",
                prg->hard_fn->n, prg->hard_fn->hardness, prg->hard_fn->epsilon);
        if (prg->hard_fn->description)
            fprintf(fp, "  description: %s\n", prg->hard_fn->description);
    }
    fprintf(fp, "  security: %zu bits\n", nw_prg_security_level(prg));
    fprintf(fp, "}\n");
}

NWPRG *nw_prg_clone(const NWPRG *prg) {
    if (!prg) return NULL;

    /* Deep copy of design */
    NWDesign *design_copy = nw_design_create(
        prg->design->k, prg->design->m, prg->design->l,
        prg->design->intersect_bound
    );
    if (!design_copy) return NULL;

    for (size_t i = 0; i < prg->design->m; i++)
        memcpy(design_copy->sets[i], prg->design->sets[i],
               prg->design->l * sizeof(size_t));

    /* Shallow copy of hard function (shared) */
    HardFunction *hf_copy = nw_hard_function_create(
        prg->hard_fn->n, prg->hard_fn->evaluate,
        prg->hard_fn->hardness, prg->hard_fn->epsilon,
        prg->hard_fn->description
    );
    if (!hf_copy) { nw_design_free(design_copy); return NULL; }

    NWPRG *clone = nw_prg_construct(design_copy, hf_copy);
    if (!clone) {
        nw_hard_function_free(hf_copy);
        nw_design_free(design_copy);
        return NULL;
    }

    return clone;
}
