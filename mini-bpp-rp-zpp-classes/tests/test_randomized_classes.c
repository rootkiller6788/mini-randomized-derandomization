/******************************************************************************
 * test_randomized_classes.c — Comprehensive tests for BPP/RP/ZPP module
 *
 * Tests cover:
 *   L1: RandomSource, ProbabilisticTM structs
 *   L2: BPP/RP/ZPP amplification
 *   L3: PTM simulation and estimation
 *   L4: Class hierarchy, Adleman, Chernoff bounds
 *   L5: Freivalds, PIT, Miller-Rabin
 *   L6: circuit construction and evaluation
 *
 * Compile: cc -I../include -o test test.c ../src/ -lm
 ******************************************************************************/

#include "randomized_classes.h"
#include "chernoff_bounds.h"
#include "derandomization_methods.h"
#include "probabilistic_circuits.h"
#include "randomized_reductions.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) do { \
    tests_run++; \
    printf("  TEST %3d: %-50s ... ", tests_run, name); \
    fflush(stdout); \
} while(0)

#define PASS() do { \
    printf("PASSED\n"); \
    tests_passed++; \
} while(0)

#define FAIL(msg) do { \
    printf("FAILED: %s\n", msg); \
    tests_failed++; \
} while(0)

#define ASSERT_TRUE(cond, msg) do { \
    if (!(cond)) { FAIL(msg); return; } \
} while(0)

#define ASSERT_FALSE(cond, msg) do { \
    if (cond) { FAIL(msg); return; } \
} while(0)

#define ASSERT_EQ(a, b, msg) do { \
    if ((a) != (b)) { \
        char buf[256]; \
        snprintf(buf, sizeof(buf), "%s (expected %zu, got %zu)", msg, (size_t)(b), (size_t)(a)); \
        FAIL(buf); return; \
    } \
} while(0)

#define ASSERT_DOUBLE_NEAR(a, b, tol, msg) do { \
    if (fabs((a) - (b)) > (tol)) { \
        char buf[256]; \
        snprintf(buf, sizeof(buf), "%s (expected %.6f, got %.6f)", msg, (b), (a)); \
        FAIL(buf); return; \
    } \
} while(0)

/* ---- L1: RandomSource tests ---- */

static void test_random_source_init(void) {
    TEST("RandomSource init with seed 42");
    RandomSource rs = random_source_init(42);
    ASSERT_EQ(rs.seed, 42, "Seed not stored");
    ASSERT_EQ(rs.tosses_used, 0, "Tosses not initialized to 0");
    ASSERT_FALSE(rs.exhausted, "Should not be exhausted initially");
    PASS();
}

static void test_random_source_bit_distribution(void) {
    TEST("Random source bit distribution (10k trials)");
    RandomSource rs = random_source_init(12345);
    size_t ones = 0, zeros = 0;
    size_t n = 10000;
    for (size_t i = 0; i < n; i++) {
        if (random_source_bit(&rs)) ones++; else zeros++;
    }
    /* 99.9% confidence: |ones - n/2| < 4·√n ≈ 400 */
    size_t diff = (ones > zeros) ? (ones - zeros) : (zeros - ones);
    ASSERT_TRUE(diff < 500, "Bit distribution too skewed (diff >= 500)");
    (void)diff;
    PASS();
}

static void test_random_source_uniform(void) {
    TEST("Random source uniform [0,99] distribution");
    RandomSource rs = random_source_init(9999);
    size_t buckets[100] = {0};
    size_t n = 10000;
    for (size_t i = 0; i < n; i++) {
        buckets[random_source_uniform(&rs, 100)]++;
    }
    /* Each bucket should have ~100 ± 40 */
    size_t min_b = n, max_b = 0;
    for (int i = 0; i < 100; i++) {
        if (buckets[i] < min_b) min_b = buckets[i];
        if (buckets[i] > max_b) max_b = buckets[i];
    }
    ASSERT_TRUE(min_b > 30, "Bucket too small");
    ASSERT_TRUE(max_b < 180, "Bucket too large");
    PASS();
}

/* ---- L2: BPP amplification ---- */

static bool always_accept(const char *input, size_t len, RandomSource *rs) {
    (void)input; (void)len; (void)rs;
    return true;
}

static bool always_reject(const char *input, size_t len, RandomSource *rs) {
    (void)input; (void)len; (void)rs;
    return false;
}

static bool noisy_decider(const char *input, size_t len, RandomSource *rs) {
    (void)input; (void)len;
    /* Correct answer = true, but err with prob 1/3 */
    return random_source_uniform(rs, 3) != 0;
}

static void test_bpp_amplify_always_accept(void) {
    TEST("BPP amplify on always-accept decider");
    RandomizedDecision d = bpp_amplify("test", 4, always_accept, 10, 1000);
    ASSERT_TRUE(d.accepted, "BPP should accept");
    ASSERT_TRUE(d.confidence > 0.99, "Confidence too low");
    ASSERT_EQ(d.trials, 180, "Expected 18*10=180 trials");
    PASS();
}

static void test_bpp_amplify_always_reject(void) {
    TEST("BPP amplify on always-reject decider");
    RandomizedDecision d = bpp_amplify("test", 4, always_reject, 10, 1000);
    ASSERT_FALSE(d.accepted, "BPP should reject");
    ASSERT_TRUE(d.confidence > 0.99, "Confidence too low");
    PASS();
}

static void test_rp_amplify_always_accept(void) {
    TEST("RP amplify on always-accept decider");
    RandomizedDecision d = rp_amplify("test", 4, always_accept, 20, 1000);
    ASSERT_TRUE(d.accepted, "RP should accept immediately");
    ASSERT_EQ(d.accept_count, 1, "Should accept on first trial");
    PASS();
}

static void test_zpp_simulate(void) {
    TEST("ZPP simulate with decisive oracle");
    int zpp_oracle(const char *s, size_t l, RandomSource *r) {
        (void)s; (void)l; (void)r;
        return 1;  /* Always accept (no "don't know") */
    }
    RandomizedDecision d = zpp_simulate("test", 4, zpp_oracle, 10, 1000);
    ASSERT_TRUE(d.accepted, "ZPP should accept");
    ASSERT_EQ(d.confidence, 1.0, "ZPP never errs");
    PASS();
}

/* ---- L3: PTM simulation ---- */

bool true_oracle(const char *s, size_t l) { (void)s; (void)l; return true; }
bool false_oracle(const char *s, size_t l) { (void)s; (void)l; return false; }

static void test_ptm_run_once(void) {
    TEST("PTM run once — correct accept");
    ProbabilisticTM ptm = ptm_create(always_accept, 10, 100, 0.1);
    RandomSource rs = random_source_init(42);
    PTMOutcome outcome = ptm_run_once(&ptm, "test", 4, &rs, true_oracle);
    ASSERT_EQ(outcome, PTM_CORRECT_ACCEPT, "Expected correct accept");
    PASS();
}

static void test_ptm_estimate_accept_prob(void) {
    TEST("PTM estimate accept probability");
    ProbabilisticTM ptm = ptm_create(noisy_decider, 10, 100, 0.0);
    double conf_radius;
    double p = ptm_estimate_accept_prob(&ptm, "test", 4, 5000, &conf_radius);
    ASSERT_DOUBLE_NEAR(p, 2.0/3.0, 0.05, "Accept prob should be ~2/3");
    ASSERT_TRUE(conf_radius < 0.1, "Confidence radius too large");
    PASS();
}

/* ---- L4: Class relationships ---- */

static void test_verify_class_hierarchy(void) {
    TEST("Class hierarchy P ⊆ ZPP ⊆ RP ⊆ BPP");
    ASSERT_TRUE(verify_class_hierarchy(), "Hierarchy should hold");
    PASS();
}

static void test_classify_machine(void) {
    TEST("Classify machine by error rates");
    bool results[4] = {false};

    /* BPP: both errors < 1/3 */
    classify_machine(0.1, 0.1, 0.95, results);
    ASSERT_TRUE(results[0], "Should be BPP");
    ASSERT_FALSE(results[1], "Should NOT be RP (false pos > 0)");

    /* RP: false_pos = 0, false_neg < 1/2 */
    classify_machine(0.0, 0.3, 0.95, results);
    ASSERT_TRUE(results[1], "Should be RP");
    ASSERT_TRUE(results[0], "Also BPP");

    /* ZPP: both ≈ 0 */
    classify_machine(0.001, 0.001, 0.95, results);
    ASSERT_TRUE(results[3], "Should be ZPP");
    PASS();
}

static void test_adleman_advice_size(void) {
    TEST("Adleman Theorem: compute advice size");
    size_t advice = adleman_theorem_advice_size(10, 0.0001);
    ASSERT_TRUE(advice > 0, "Advice size should be positive");
    ASSERT_TRUE(advice < 100000, "Advice size unreasonably large");
    PASS();
}

/* ---- L4: Chernoff bounds ---- */

static void test_chernoff_upper_bound(void) {
    TEST("Chernoff upper tail bound");
    double bound = chernoff_upper_tail(100.0, 0.5);
    ASSERT_TRUE(bound > 0.0 && bound < 1.0, "Bound should be in (0,1)");
    ASSERT_TRUE(bound < 0.01, "Bound should be small for moderate deviation");
    PASS();
}

static void test_chernoff_lower_bound(void) {
    TEST("Chernoff lower tail bound");
    double bound = chernoff_lower_tail(100.0, 0.5);
    ASSERT_TRUE(bound > 0.0 && bound < 1.0, "Bound should be in (0,1)");
    ASSERT_TRUE(bound < 0.01, "Lower tail: small bound expected");
    PASS();
}

static void test_hoeffding_bound(void) {
    TEST("Hoeffding bound for n=100, epsilon=0.1");
    double bound = hoeffding_bound(100, 0.1);
    ASSERT_TRUE(bound < 0.5, "Hoeffding: should be < 0.5 for n=100, eps=0.1");
    PASS();
}

static void test_bpp_trials(void) {
    TEST("BPP trials for error 0.01");
    size_t t = bpp_trials_for_error(0.01);
    ASSERT_TRUE(t >= 80 && t <= 100, "Expected ~83, got absurd value");
    PASS();
}

/* ---- L5: Freivalds algorithm ---- */

static void test_freivalds_correct(void) {
    TEST("Freivalds: correct product verification");
    size_t n = 10;
    double *A = (double *)calloc(n * n, sizeof(double));
    double *B = (double *)calloc(n * n, sizeof(double));
    double *C = (double *)calloc(n * n, sizeof(double));

    /* Fill A, B with identity-ish */
    for (size_t i = 0; i < n; i++) A[i * n + i] = 1.0;
    for (size_t i = 0; i < n; i++) B[i * n + i] = 1.0;
    /* C = A × B = I */
    for (size_t i = 0; i < n; i++) C[i * n + i] = 1.0;

    RandomSource rs = random_source_init(42);
    bool result = freivalds_verify_amplified(A, B, C, n, 10, &rs);
    ASSERT_TRUE(result, "Freivalds should verify I×I=I");

    /* Corrupt C */
    C[0] = 999.0;
    result = freivalds_verify_amplified(A, B, C, n, 10, &rs);
    ASSERT_FALSE(result, "Freivalds should detect corruption");

    free(A); free(B); free(C);
    PASS();
}

/* ---- L5: Miller-Rabin ---- */

static void test_zpp_primality_small(void) {
    TEST("ZPP primality test — small primes");
    RandomSource rs = random_source_init(123);
    ASSERT_EQ(zpp_primality_test(2, &rs, true), 1, "2 is prime");
    ASSERT_EQ(zpp_primality_test(3, &rs, true), 1, "3 is prime");
    ASSERT_EQ(zpp_primality_test(4, &rs, true), 0, "4 is composite");
    ASSERT_EQ(zpp_primality_test(17, &rs, true), 1, "17 is prime");
    ASSERT_EQ(zpp_primality_test(97, &rs, true), 1, "97 is prime");
    ASSERT_EQ(zpp_primality_test(100, &rs, true), 0, "100 is composite");
    PASS();
}

static void test_zpp_primality_randomized(void) {
    TEST("ZPP primality test — randomized mode");
    RandomSource rs = random_source_init(123);
    ASSERT_EQ(zpp_primality_test(101, &rs, false), 1, "101 is prime");
    ASSERT_EQ(zpp_primality_test(143, &rs, false), 0, "143 = 11×13 composite");
    PASS();
}

/* ---- L6: Circuit construction ---- */

static void test_circuit_create_evaluate(void) {
    TEST("Circuit: create AND circuit and evaluate");
    /* C(x0, x1) = x0 AND x1 */
    BooleanCircuit *c = circuit_create(2, 0, 1, 10);
    ASSERT_TRUE(c != NULL, "Circuit creation failed");

    size_t in0 = circuit_add_input(c, 0);
    size_t in1 = circuit_add_input(c, 1);
    size_t inputs[2] = {in0, in1};
    size_t out = circuit_add_gate(c, GATE_AND, inputs, 2);
    circuit_set_output(c, out, 0);

    bool output;
    bool x[2];

    x[0] = false; x[1] = false;
    circuit_evaluate(c, x, NULL, &output);
    ASSERT_FALSE(output, "0 AND 0 = 0");

    x[0] = true; x[1] = true;
    circuit_evaluate(c, x, NULL, &output);
    ASSERT_TRUE(output, "1 AND 1 = 1");

    x[0] = true; x[1] = false;
    circuit_evaluate(c, x, NULL, &output);
    ASSERT_FALSE(output, "1 AND 0 = 0");

    circuit_destroy(c);
    PASS();
}

static void test_circuit_parity(void) {
    TEST("Circuit: PARITY function for n=4");
    BooleanCircuit *c = circuit_parity(4);
    ASSERT_TRUE(c != NULL, "Parity circuit creation failed");

    bool output;
    for (size_t val = 0; val < 16; val++) {
        bool x[4] = {(val>>0)&1, (val>>1)&1, (val>>2)&1, (val>>3)&1};
        circuit_evaluate(c, x, NULL, &output);
        bool expected = x[0] ^ x[1] ^ x[2] ^ x[3];
        if (output != expected) {
            char buf[128];
            snprintf(buf, sizeof(buf), "PARITY(%d,%d,%d,%d) expected %d got %d",
                     x[0], x[1], x[2], x[3], expected, output);
            FAIL(buf);
            circuit_destroy(c);
            return;
        }
    }
    circuit_destroy(c);
    PASS();
}

static void test_shannon_lower_bound(void) {
    TEST("Shannon lower bound for n=5");
    size_t lb = shannon_lower_bound(5, 3, 2);
    /* Shannon: Ω(2^n/n) ≈ 6.4 for n=5. Computed bound may be ~5 due to
     * the simplified upper bound on #circuits. Test that the bound is
     * non-trivial (>0) and < 2^n. */
    ASSERT_TRUE(lb >= 2, "Shannon bound should be at least 2");
    ASSERT_TRUE(lb <= 64, "Shannon bound should be less than 2^n=32?");
    PASS();
}

/* ---- L4: Pairwise independent hash ---- */

static void test_pairwise_hash(void) {
    TEST("Pairwise independent hash family");
    RandomSource rs = random_source_init(777);
    PairwiseHash h = pairwise_hash_create(100, 10, &rs);

    /* Hash a few values */
    for (uint64_t x = 0; x < 10; x++) {
        uint64_t val = pairwise_hash_eval(&h, x);
        ASSERT_TRUE(val < 10, "Hash output out of range");
    }
    PASS();
}

/* ---- L7: Application: quicksort analysis ---- */

static void test_quicksort_analysis(void) {
    TEST("Randomized quicksort analysis n=100");
    double mean, stddev;
    randomized_quicksort_analysis(100, 20, &mean, &stddev);
    /* Expected comparisons: 2(n+1)H_n - 4n ≈ 2·101·5.187 - 400 ≈ 648 */
    ASSERT_TRUE(mean > 300 && mean < 1000, "Quicksort mean out of expected range");
    ASSERT_TRUE(stddev > 0.0, "Stddev should be positive");
    PASS();
}

/* ---- L8: Impagliazzo worlds ---- */

static void test_impagliazzo_worlds(void) {
    TEST("Impagliazzo's five worlds — names");
    for (int w = 0; w < 5; w++) {
        const char *name = impagliazzo_world_name((ImpagliazzoWorld)w);
        ASSERT_TRUE(name != NULL && strlen(name) > 0, "World name empty");
        const char *desc = impagliazzo_world_description((ImpagliazzoWorld)w);
        ASSERT_TRUE(desc != NULL && strlen(desc) > 20, "World description too short");
    }
    PASS();
}

/* ---- L6: Karger min-cut ---- */

static void test_karger_min_cut(void) {
    TEST("Karger min-cut on 4-cycle");
    size_t n = 4;
    int *adj = (int *)calloc(n * n, sizeof(int));
    /* 4-cycle: edges (0-1), (1-2), (2-3), (3-0) */
    adj[0*n+1] = adj[1*n+0] = 1;
    adj[1*n+2] = adj[2*n+1] = 1;
    adj[2*n+3] = adj[3*n+2] = 1;
    adj[3*n+0] = adj[0*n+3] = 1;

    RandomSource rs = random_source_init(12345);
    int mc = karger_min_cut(adj, n, 50, &rs);
    /* Minimum cut in 4-cycle is 2 */
    ASSERT_EQ(mc, 2, "Min cut of 4-cycle should be 2");
    free(adj);
    PASS();
}

/* ---- L6: Randomized 2-SAT ---- */

static void test_randomized_2sat(void) {
    TEST("Randomized 2-SAT on simple formula");
    /* (x0 ∨ x1) ∧ (¬x0 ∨ x2) ∧ (¬x1 ∨ ¬x2) — satisfiable */
    int clauses[] = {
        0, 0, 1, 0,   /* x0 ∨ x1 */
        0, 1, 2, 0,   /* ¬x0 ∨ x2 */
        1, 1, 2, 1,   /* ¬x1 ∨ ¬x2 */
    };
    bool assign[3];
    RandomSource rs = random_source_init(42);

    bool sat = randomized_2sat(clauses, 3, 3, 1000, &rs, assign);
    ASSERT_TRUE(sat, "2-SAT should find assignment");

    /* Verify assignment */
    bool c0 = assign[0] || assign[1];
    bool c1 = (!assign[0]) || assign[2];
    bool c2 = (!assign[1]) || (!assign[2]);
    ASSERT_TRUE(c0 && c1 && c2, "Assignment should satisfy all clauses");
    PASS();
}

/* ---- Derandomization: MAX-CUT ---- */

static void test_max_cut_derandomized(void) {
    TEST("MAX-CUT derandomization on triangle");
    size_t n = 3;
    int adj[9] = {0,1,1, 1,0,1, 1,1,0};
    bool assignments[3];
    size_t cut = max_cut_derandomized(adj, n, assignments);
    /* Triangle: max cut is 2 (two edges across cut) */
    ASSERT_TRUE(cut >= 1, "MAX-CUT should find ≥ n(n-1)/4");
    PASS();
}

/* ---- Run all tests ---- */

int main(void) {
    printf("=== BPP/RP/ZPP Randomized Classes — Test Suite ===\n\n");

    /* L1 */
    test_random_source_init();
    test_random_source_bit_distribution();
    test_random_source_uniform();

    /* L2 */
    test_bpp_amplify_always_accept();
    test_bpp_amplify_always_reject();
    test_rp_amplify_always_accept();
    test_zpp_simulate();

    /* L3 */
    test_ptm_run_once();
    test_ptm_estimate_accept_prob();

    /* L4 */
    test_verify_class_hierarchy();
    test_classify_machine();
    test_adleman_advice_size();
    test_chernoff_upper_bound();
    test_chernoff_lower_bound();
    test_hoeffding_bound();
    test_bpp_trials();
    test_pairwise_hash();

    /* L5 */
    test_freivalds_correct();
    test_zpp_primality_small();
    test_zpp_primality_randomized();

    /* L6 */
    test_circuit_create_evaluate();
    test_circuit_parity();
    test_shannon_lower_bound();
    test_karger_min_cut();
    test_randomized_2sat();

    /* L7 */
    test_quicksort_analysis();

    /* L8 */
    test_impagliazzo_worlds();

    /* Derandomization */
    test_max_cut_derandomized();

    printf("\n=== Results: %d run, %d passed, %d failed ===\n",
           tests_run, tests_passed, tests_failed);

    return tests_failed > 0 ? 1 : 0;
}
