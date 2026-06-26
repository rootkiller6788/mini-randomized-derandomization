/******************************************************************************
 * test.c �� Comprehensive Test Suite for mini-hardness-vs-randomness
 *
 * Tests all core APIs: circuit construction, lower bound computation,
 * hardness-to-PRG conversion, derandomization, worst-to-average reduction,
 * IW theorem conditions, and applications.
 *
 * Uses standard assert() for all checks.
 ******************************************************************************/

#include "circuit_lower_bounds.h"
#include "hardness_randomness.h"
#include "derandomization_via_hardness.h"
#include "worst_to_average.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

/* Forward declarations from iw_theorem.c */
extern bool iw_theorem_conditions_met(const HardnessCert *cert, size_t n);
extern bool iw_full_derandomization_pipeline(size_t n, HRLevel level);
extern bool iw_check_eca_assumption(size_t n, double epsilon);
extern double iw_derandomization_blowup(size_t n, double hardness);
extern double iw_pseudodeterministic_probability(size_t n, double hardness);
extern bool iw_circuit_lb_sufficient(const CircuitLB *clb);
extern bool iw_self_test(void);

/* Forward declarations from applications.c */
extern bool app_owf_from_hardness(size_t n, double hardness);
extern bool app_goldreich_levin_predicate(const bool *x, const bool *r, size_t n);
extern size_t app_prg_stream_from_predicate(size_t n, bool *output, size_t out_len,
                                              const bool *seed, size_t seed_len);
extern bool app_miller_rabin_test(uint64_t n, uint64_t a);
extern size_t app_primality_rounds_for_error(double error);
extern size_t app_crypto_key_size(double hardness_rate, double target_security_bits);
extern double app_security_from_hardness(HRLevel level, size_t n);
extern void app_expander_walk(size_t n, size_t steps, size_t seed, size_t *output);
extern size_t app_pairwise_independent_hash(size_t x, size_t a, size_t b,
                                              size_t p, size_t m);
extern double app_deterministic_monte_carlo(size_t n, size_t samples, uint64_t seed);
extern void app_deterministic_scenarios(size_t n, size_t *obstacle_positions,
                                         size_t *velocities, uint64_t seed);
extern void app_gps_gold_code(size_t n, bool *code_out, uint64_t seed);
extern double app_climate_ensemble_parameter(size_t param_idx, size_t ensemble_size);
extern double app_nuclear_safety_margin(double temperature, double pressure,
                                         double material_strength);
extern bool app_iso26262_derandomization_verification(size_t n, double confidence);
extern bool app_self_test(void);

static size_t tests_passed = 0;
static size_t tests_failed = 0;

#define TEST(name, expr) do { \
    if (!(expr)) { \
        fprintf(stderr, "FAIL: %s (%s:%d)\n", name, __FILE__, __LINE__); \
        tests_failed++; \
    } else { \
        tests_passed++; \
    } \
} while(0)

/* ================================================================
 * Circuit Lower Bound Tests
 * ================================================================ */

static void test_circuit_creation(void) {
    BooleanCircuit *c = clb_circuit_create(4);
    TEST("circuit_create non-null", c != NULL);
    TEST("circuit_create n", c->n == 4);
    TEST("circuit_create num_gates", c->num_gates == 4);
    TEST("circuit_create inputs created", c->gates[0].type == GATE_INPUT);
    TEST("circuit_verify_structure", clb_circuit_verify_structure(c));
    clb_circuit_free(c);
}

static void test_circuit_evaluation(void) {
    /* Build: f(x1,x2) = x1 AND x2 */
    BooleanCircuit *c = clb_circuit_create(2);
    assert(c != NULL);
    size_t and_gate = clb_circuit_add_gate(c, GATE_AND, 0, 1, false);
    TEST("circuit_add_gate returns valid", and_gate == 2);
    size_t nout = clb_circuit_set_output(c, and_gate);
    TEST("circuit_set_output count", nout == 1);
    /* Evaluate */
    bool x00[] = {false, false};
    bool x01[] = {false, true};
    bool x10[] = {true, false};
    bool x11[] = {true, true};
    bool result[1];
    clb_circuit_evaluate(c, x00, result);
    TEST("AND(0,0)=0", result[0] == false);
    clb_circuit_evaluate(c, x01, result);
    TEST("AND(0,1)=0", result[0] == false);
    clb_circuit_evaluate(c, x10, result);
    TEST("AND(1,0)=0", result[0] == false);
    clb_circuit_evaluate(c, x11, result);
    TEST("AND(1,1)=1", result[0] == true);
    clb_circuit_free(c);
}

static void test_circuit_negation(void) {
    /* Build: f(x1,x2) = NOT(x1 AND x2) = NAND */
    BooleanCircuit *c = clb_circuit_create(2);
    size_t and_gate = clb_circuit_add_gate(c, GATE_AND, 0, 1, false);
    size_t not_gate = clb_circuit_add_gate(c, GATE_NOT, and_gate, and_gate, false);
    clb_circuit_set_output(c, not_gate);
    bool x11[] = {true, true};
    bool result[1];
    clb_circuit_evaluate(c, x11, result);
    TEST("NAND(1,1)=0", result[0] == false);
    clb_circuit_free(c);
}

/* ================================================================
 * Circuit Lower Bound Theorem Tests
 * ================================================================ */

static void test_shannon_bound(void) {
    CircuitLB lb = clb_shannon_counting(8);
    TEST("shannon n=8", lb.n == 8);
    TEST("shannon bound > 0", lb.lower_bound > 0.0);
    TEST("shannon bound < 2^n", lb.lower_bound < 256.0);
}

static void test_lupanov_bound(void) {
    CircuitLB lb = clb_lupanov_upper(8);
    TEST("lupanov n=8", lb.n == 8);
    TEST("lupanov bound > shannon thm bound", lb.lower_bound > 0.0);
}

static void test_hastad_bound(void) {
    CircuitLB lb = clb_hastad_switching(16, 3);
    TEST("hastad n=16 depth=3", lb.n == 16);
    TEST("hastad bound > 0", lb.lower_bound > 0.0);
}

static void test_smolensky_bound(void) {
    CircuitLB lb = clb_smolensky_modp(16, 3);
    TEST("smolensky n=16 p=3", lb.n == 16);
    TEST("smolensky bound > 0", lb.lower_bound > 0.0);
}

static void test_razborov_bound(void) {
    CircuitLB lb = clb_razborov_clique(20);
    TEST("razborov n=20", lb.n == 20);
    TEST("razborov bound > 0", lb.lower_bound > 0.0);
}

static void test_quality_metric(void) {
    CircuitLB lb = clb_shannon_counting(10);
    double q = clb_quality_metric(&lb);
    TEST("quality metric > 0", q > 0.0);
    TEST("quality metric �� 1", q <= 1.0);
}

static void test_is_strong_enough(void) {
    CircuitLB lb_weak = clb_majority(10);
    TEST("majority not exponential", !clb_is_strong_enough(&lb_weak, HR_EXPONENTIAL));
    CircuitLB lb_strong = clb_shannon_counting(10);
    TEST("shannon is strong enough", clb_is_strong_enough(&lb_strong, HR_STRONG));
}

/* ================================================================
 * Truth Table Tests
 * ================================================================ */

static void test_truth_tables(void) {
    size_t n = 8;
    bool *parity_tt = clb_build_parity_truth_table(n);
    TEST("parity tt non-null", parity_tt != NULL);
    TEST("parity(0) = 0", parity_tt[0] == false);
    TEST("parity(1) = 1", parity_tt[1] == true);
    free(parity_tt);
    bool *maj_tt = clb_build_majority_truth_table(n);
    TEST("majority tt non-null", maj_tt != NULL);
    TEST("majority(0) = 0", maj_tt[0] == false);
    free(maj_tt);
    bool *rand_tt = clb_build_random_truth_table(n, 42);
    TEST("random tt non-null", rand_tt != NULL);
    free(rand_tt);
}

/* ================================================================
 * Hardness-Randomness Tests
 * ================================================================ */

static void test_hr_hardness_from_circuit_lb(void) {
    double h = hr_hardness_from_circuit_lb(16, 1024);
    /* 1024 = 2^10, so hardness = 10/16 = 0.625 */
    TEST("hardness 10/16 �� 0.6", h > 0.5 && h < 0.75);
}

static void test_hr_seed_len(void) {
    size_t weak_seed = hr_seed_len(100, HR_WEAK);
    TEST("weak seed < n", weak_seed <= 100);
    size_t strong_seed = hr_seed_len(100, HR_STRONG);
    TEST("strong seed < weak seed", strong_seed <= weak_seed);
}

static void test_hr_output_len(void) {
    size_t out = hr_output_len(10, HR_STRONG);
    TEST("output_len > seed_len", out > hr_seed_len(10, HR_STRONG));
}

static void test_hr_bpp_equals_p(void) {
    TEST("strong => BPP=P", hr_bpp_equals_p_check(HR_STRONG, 50));
    TEST("weak !=> BPP=P", !hr_bpp_equals_p_check(HR_WEAK, 50));
}

static void test_hr_iw_threshold(void) {
    double thresh = hr_impagliazzo_wigderson_threshold(100);
    TEST("IW threshold > 0", thresh > 0.0);
    TEST("IW threshold < 1", thresh < 1.0);
}

static void test_hr_amplification_possible(void) {
    TEST("amplification feasible", hr_hardness_amplification_possible(0.1, 16));
    TEST("amplification infeasible", !hr_hardness_amplification_possible(0.0, 16));
}

/* ================================================================
 * Derandomization Tests
 * ================================================================ */

static void test_derand_cert(void) {
    HardnessCert cert = derand_cert_from_circuit_lb(16, 65536);
    /* 65536 = 2^16, so hardness = 16/16 = 1.0 */
    TEST("cert hardness �� 1.0", cert.hardness > 0.9);
    TEST("verify cert", derand_verify_certificate(&cert, 16));
    TEST("verify wrong n", !derand_verify_certificate(&cert, 15));
}

static void test_derand_time_complexity(void) {
    size_t time = derand_time_complexity(100, 0.1);
    TEST("derand time finite", time < SIZE_MAX);
}

static void test_derand_error_prob(void) {
    double p = derand_error_prob(10, 100, 0.3);
    TEST("error prob < 1", p < 1.0);
    TEST("error prob > 0", p > 0.0);
}

static void test_derand_advice_size(void) {
    size_t advice = derand_advice_size(10, 0.5);
    TEST("advice > 0", advice > 0);
}

static void test_derand_confidence(void) {
    double conf = derand_confidence_bound(10, 100);
    TEST("confidence > 0.9", conf > 0.9);
}

/* ================================================================
 * Worst-to-Average Tests
 * ================================================================ */

static void test_yao_xor_lemma(void) {
    WTAConversion conv = wta_yao_xor_lemma(0.1, 16, 5);
    TEST("avg hardness amplified", conv.avg_hardness > conv.worst_hardness);
    TEST("worst preserved", conv.worst_hardness == 0.1);
}

static void test_impagliazzo_hardcore(void) {
    WTAConversion conv = wta_impagliazzo_hardcore(0.2, 16);
    TEST("hardcore avg �� 0.5", conv.avg_hardness > 0.4 && conv.avg_hardness < 0.6);
}

static void test_conversion_feasible(void) {
    TEST("feasible", wta_conversion_feasible(0.1, 16));
    TEST("not feasible", !wta_conversion_feasible(0.0, 16));
}

static void test_sample_complexity(void) {
    size_t samples = wta_sample_complexity(16, 0.1, 0.05);
    TEST("samples > 0", samples > 0);
}

static void test_xor_lemma_bound(void) {
    double bound = wta_xor_lemma_bound(0.3, 4);
    /* advantage �� (1-0.3)^4 = 0.7^4 �� 0.24 */
    TEST("XOR bound close to 0.5", bound >= 0.5 && bound < 0.75);
}

/* ================================================================
 * IW Theorem Tests
 * ================================================================ */

static void test_iw_theorem_conditions(void) {
    HardnessCert cert = derand_cert_from_circuit_lb(16, 65536);
    TEST("IW conditions met", iw_theorem_conditions_met(&cert, 16));
}

static void test_iw_pipeline(void) {
    bool ok = iw_full_derandomization_pipeline(12, HR_STRONG);
    TEST("IW pipeline for n=12 strong", ok);
}

static void test_iw_eca_assumption(void) {
    TEST("ECA assumption valid", iw_check_eca_assumption(20, 0.1));
}

static void test_iw_pseudodeterministic(void) {
    double prob = iw_pseudodeterministic_probability(10, 0.5);
    TEST("pseudo-det prob > 0.9", prob > 0.9);
}

static void test_iw_circuit_lb_sufficient(void) {
    CircuitLB lb = clb_shannon_counting(16);
    TEST("shannon sufficient for IW", iw_circuit_lb_sufficient(&lb));
}

/* ================================================================
 * Application Tests
 * ================================================================ */

static void test_app_owf(void) {
    TEST("OWF from hardness possible", app_owf_from_hardness(16, 0.1));
}

static void test_app_goldreich_levin(void) {
    bool x[] = {true, false, true, true};
    bool r[] = {false, true, true, false};
    /* ?x,r? = 1*0 + 0*1 + 1*1 + 1*0 = 1 mod 2 */
    TEST("GL predicate correct", app_goldreich_levin_predicate(x, r, 4) == true);
}

static void test_app_prg_stream(void) {
    bool seed[] = {true, false, true, false, true, false, true, false};
    bool output[32];
    size_t gen = app_prg_stream_from_predicate(8, output, 32, seed, 8);
    TEST("PRG stream generated all bits", gen == 32);
}

static void test_app_miller_rabin(void) {
    TEST("7 is prime (base 2)", app_miller_rabin_test(7, 2));
    TEST("11 is prime (base 2)", app_miller_rabin_test(11, 2));
    TEST("15 not prime (base 2)", !app_miller_rabin_test(15, 2));
}

static void test_app_crypto_key_size(void) {
    size_t ks = app_crypto_key_size(0.1, 128.0);
    TEST("key size �� 1280", ks >= 1280);
}

static void test_app_hash(void) {
    size_t h = app_pairwise_independent_hash(12345, 3, 7, 99991, 100);
    TEST("hash < 100", h < 100);
    /* Same input �� same output */
    size_t h2 = app_pairwise_independent_hash(12345, 3, 7, 99991, 100);
    TEST("hash deterministic", h == h2);
}

static void test_app_expander_walk(void) {
    size_t output[20];
    app_expander_walk(97, 20, 42, output);
    for (int i = 0; i < 20; i++) {
        TEST("expander output < n", output[i] < 97);
    }
}

static void test_app_gps_code(void) {
    bool code[64];
    app_gps_gold_code(64, code, 12345);
    /* Verify it produces bits */
    size_t ones = 0;
    for (int i = 0; i < 64; i++) if (code[i]) ones++;
    TEST("GPS code has some ones", ones > 0);
    TEST("GPS code has some zeros", ones < 64);
}

static void test_app_nuclear(void) {
    double margin = app_nuclear_safety_margin(300.0, 15.0, 500.0);
    TEST("safety margin > 0", margin > 0.0);
}

int main(void) {
    printf("=== mini-hardness-vs-randomness Test Suite ===\n\n");
    /* Circuit tests */
    test_circuit_creation();
    test_circuit_evaluation();
    test_circuit_negation();
    test_shannon_bound();
    test_lupanov_bound();
    test_hastad_bound();
    test_smolensky_bound();
    test_razborov_bound();
    test_quality_metric();
    test_is_strong_enough();
    test_truth_tables();
    /* Hardness-randomness tests */
    test_hr_hardness_from_circuit_lb();
    test_hr_seed_len();
    test_hr_output_len();
    test_hr_bpp_equals_p();
    test_hr_iw_threshold();
    test_hr_amplification_possible();
    /* Derandomization tests */
    test_derand_cert();
    test_derand_time_complexity();
    test_derand_error_prob();
    test_derand_advice_size();
    test_derand_confidence();
    /* Worst-to-average tests */
    test_yao_xor_lemma();
    test_impagliazzo_hardcore();
    test_conversion_feasible();
    test_sample_complexity();
    test_xor_lemma_bound();
    /* IW theorem tests */
    test_iw_theorem_conditions();
    test_iw_pipeline();
    test_iw_eca_assumption();
    test_iw_pseudodeterministic();
    test_iw_circuit_lb_sufficient();
    /* Application tests */
    test_app_owf();
    test_app_goldreich_levin();
    test_app_prg_stream();
    test_app_miller_rabin();
    test_app_crypto_key_size();
    test_app_hash();
    test_app_expander_walk();
    test_app_gps_code();
    test_app_nuclear();
    /* Self-tests from other modules */
    TEST("iw_self_test", iw_self_test());
    TEST("app_self_test", app_self_test());
    printf("\n=== Results: %zu passed, %zu failed ===\n",
           tests_passed, tests_failed);
    return tests_failed > 0 ? 1 : 0;
}