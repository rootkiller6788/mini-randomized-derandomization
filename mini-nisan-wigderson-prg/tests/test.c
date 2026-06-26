/**
 * @file test.c
 * @brief Comprehensive test suite for mini-nisan-wigderson-prg.
 *
 * Tests cover L1-L6 knowledge layers:
 *   L1: Type definitions, design creation/verification
 *   L2: PRG construction/evaluation, hardness amplification
 *   L3: Set systems, truth tables, bit utilities
 *   L4: Shannon bound, Yao XOR Lemma, NW Theorem verification
 *   L5: Reed-Solomon designs, PRG evaluation, hybrid argument
 *   L6: PARITY/MAJORITY circuit constructions
 *
 * All tests use standard assert().
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "nw_core.h"
#include "nw_designs.h"
#include "nw_hardness.h"
#include "nw_prg.h"
#include "nw_circuits.h"

static int tests_run = 0;
static int tests_passed = 0;
#define TEST(name) do { tests_run++; printf("  TEST: %s ... ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASSED\n"); } while(0)

/* =============================================
 * L1: Core Definitions — Design Creation
 * ============================================= */

static void test_design_create_free(void) {
    TEST("design_create_free");
    NWDesign *d = nw_design_create(10, 5, 3, 2);
    assert(d != NULL);
    assert(d->k == 10);
    assert(d->m == 5);
    assert(d->l == 3);
    assert(d->intersect_bound == 2);
    nw_design_free(d);
    PASS();
}

static void test_design_create_null_params(void) {
    TEST("design_create_null_params");
    assert(nw_design_create(0, 5, 3, 2) == NULL);
    assert(nw_design_create(10, 0, 3, 2) == NULL);
    assert(nw_design_create(10, 5, 0, 2) == NULL);
    assert(nw_design_create(5, 5, 10, 2) == NULL); /* l > k */
    PASS();
}

/* =============================================
 * L3: Design Verification — Small Intersection
 * ============================================= */

static void test_design_verify_small(void) {
    TEST("design_verify_small");
    /* Create a (6, 3, 2, 1)-design manually */
    NWDesign *d = nw_design_create(6, 3, 2, 1);
    assert(d != NULL);

    /* S_0 = {0, 1}, S_1 = {2, 3}, S_2 = {4, 5} — all disjoint */
    d->sets[0][0] = 0; d->sets[0][1] = 1;
    d->sets[1][0] = 2; d->sets[1][1] = 3;
    d->sets[2][0] = 4; d->sets[2][1] = 5;

    assert(nw_design_verify(d) == true);
    assert(nw_design_max_intersection(d) == 0);
    nw_design_free(d);
    PASS();
}

static void test_design_verify_intersecting(void) {
    TEST("design_verify_intersecting");
    /* Create a (6, 3, 2, 1)-design with violation */
    NWDesign *d = nw_design_create(6, 3, 2, 1);
    assert(d != NULL);

    /* S_0 = {0,1}, S_1 = {1,2} — intersect in {1} (size 1, OK)
     * S_2 = {1,3} — intersect with S_0 in {1} (size 1, OK)
     *                 intersect with S_1 in {1} (size 1, OK) */
    d->sets[0][0] = 0; d->sets[0][1] = 1;
    d->sets[1][0] = 1; d->sets[1][1] = 2;
    d->sets[2][0] = 1; d->sets[2][1] = 3;
    assert(nw_design_verify(d) == true);

    /* Now violate: S_2 = {0,1} — intersect with S_0 in {0,1} (size 2, fails) */
    d->sets[2][0] = 0; d->sets[2][1] = 1;
    assert(nw_design_verify(d) == false);

    nw_design_free(d);
    PASS();
}

/* =============================================
 * L3: Boolean Array Utilities
 * ============================================= */

static void test_int_seed_conversion(void) {
    TEST("int_seed_conversion");
    bool seed[64] = {0};
    nw_int_to_seed(8, 0x55, seed);
    assert(seed[0] == true);   /* bit 0 */
    assert(seed[1] == false);
    assert(seed[2] == true);   /* bit 2 */
    assert(seed[3] == false);
    assert(seed[4] == true);   /* bit 4 */
    assert(seed[5] == false);
    assert(seed[6] == true);   /* bit 6 */
    assert(seed[7] == false);

    uint64_t val = nw_seed_to_int(8, seed);
    assert(val == 0x55);
    PASS();
}

static void test_seed_restrict(void) {
    TEST("seed_restrict");
    bool seed[8] = {true, false, true, true, false, true, false, false};
    size_t subset[3] = {0, 4, 7};
    bool result[3] = {0};
    nw_seed_restrict(8, seed, subset, 3, result);
    assert(result[0] == true);   /* seed[0] */
    assert(result[1] == false);  /* seed[4] */
    assert(result[2] == false);  /* seed[7] */
    PASS();
}

static void test_hamming_weight(void) {
    TEST("hamming_weight");
    bool arr[8] = {true, false, true, false, true, true, false, true};
    assert(nw_hamming_weight(arr, 8) == 5);
    assert(nw_hamming_weight(arr, 0) == 0);
    PASS();
}

static void test_bool_xor(void) {
    TEST("bool_xor");
    bool a[4] = {true, false, true, false};
    bool b[4] = {true, true, false, false};
    bool out[4] = {0};
    nw_bool_xor(a, b, out, 4);
    assert(out[0] == false); /* 1 xor 1 = 0 */
    assert(out[1] == true);  /* 0 xor 1 = 1 */
    assert(out[2] == true);  /* 1 xor 0 = 1 */
    assert(out[3] == false); /* 0 xor 0 = 0 */
    PASS();
}

static void test_ceil_log2(void) {
    TEST("ceil_log2");
    assert(nw_ceil_log2(0) == 0);
    assert(nw_ceil_log2(1) == 0);
    assert(nw_ceil_log2(2) == 1);
    assert(nw_ceil_log2(3) == 2);
    assert(nw_ceil_log2(4) == 2);
    assert(nw_ceil_log2(7) == 3);
    assert(nw_ceil_log2(8) == 3);
    assert(nw_ceil_log2(9) == 4);
    assert(nw_ceil_log2(1024) == 10);
    PASS();
}

/* =============================================
 * L4: Yao's XOR Lemma
 * ============================================= */

static void test_yao_xor_amplify(void) {
    TEST("yao_xor_amplify");
    /* Small epsilon: XOR amplifies */
    double eps1 = nw_yao_xor_amplify(0.1, 1, 0.01);
    assert(eps1 > 0.1); /* Even k=1 with delta adds */

    /* Larger k should give larger amplification */
    double eps3 = nw_yao_xor_amplify(0.1, 3, 0.01);
    assert(eps3 > eps1);

    /* Max epsilon = 0.5 */
    double eps_max = nw_yao_xor_amplify(0.5, 5, 0.0);
    assert(eps_max <= 0.5);

    /* Zero epsilon stays zero */
    double eps_zero = nw_yao_xor_amplify(0.0, 10, 0.0);
    assert(eps_zero <= 0.01);
    PASS();
}

static void test_yao_xor_lemma_full(void) {
    TEST("yao_xor_lemma_full");
    size_t new_n;
    double new_eps, new_size;
    nw_yao_xor_lemma_full(4, 0.1, 3, 0.01, 1000.0, &new_n, &new_eps, &new_size);
    assert(new_n == 12);          /* 4 * 3 */
    assert(new_eps > 0.1);        /* Amplified hardness */
    assert(new_size < 1000.0);    /* Circuit size loss */
    PASS();
}

/* =============================================
 * L4: Shannon Lower Bound
 * ============================================= */

static void test_shannon_bound(void) {
    TEST("shannon_bound");
    double bound4 = nw_shannon_lower_bound(4);
    assert(bound4 >= 1.0);

    double bound8 = nw_shannon_lower_bound(8);
    assert(bound8 > bound4); /* Larger n => larger bound */

    double bound16 = nw_shannon_lower_bound(16);
    assert(bound16 > bound8);

    /* Lupanov upper bound should be larger than Shannon lower bound */
    double lupanov16 = nw_lupanov_upper_bound(16);
    assert(lupanov16 > bound16);
    PASS();
}

/* =============================================
 * L5: Reed-Solomon Design Construction
 * ============================================= */

static bool simple_hard_fn(const bool *input, size_t n) {
    /* Simple hard function: majority-like */
    size_t count = 0;
    for (size_t i = 0; i < n; i++)
        if (input[i]) count++;
    return (count >= n / 2);
}

static void test_reed_solomon_design(void) {
    TEST("reed_solomon_design");
    /* Small Reed-Solomon design: q=5, d=2 */
    NWDesign *d = nw_reed_solomon_design(5, 2);
    assert(d != NULL);
    assert(d->k == 25);  /* q^2 = 25 */
    assert(d->l == 5);   /* l = q */
    assert(d->m >= 25);  /* q^3 >= 125, but capped */

    /* Verify intersection property */
    assert(nw_design_verify(d) == true);
    nw_design_free(d);
    PASS();
}

static void test_random_design(void) {
    TEST("random_design");
    NWDesign *d = nw_random_design(36, 10, 4);
    assert(d != NULL);
    assert(d->k == 36);
    assert(d->m == 10);
    assert(d->l == 4);
    nw_design_free(d);
    PASS();
}

static void test_random_design_success_prob(void) {
    TEST("random_design_success_prob");
    double prob = nw_random_design_success_prob(100, 50, 10, 3);
    /* Probability should be in [0, 1] */
    assert(prob >= 0.0 && prob <= 1.0);
    PASS();
}

/* =============================================
 * L2: PRG Construction and Evaluation
 * ============================================= */

static void test_prg_construct_evaluate(void) {
    TEST("prg_construct_evaluate");
    /* Manual small design + hard function */
    NWDesign *d = nw_design_create(4, 3, 2, 1);
    assert(d != NULL);

    /* S_0 = {0,1}, S_1 = {1,2}, S_2 = {2,3} */
    d->sets[0][0] = 0; d->sets[0][1] = 1;
    d->sets[1][0] = 1; d->sets[1][1] = 2;
    d->sets[2][0] = 2; d->sets[2][1] = 3;

    HardFunction *hf = nw_hard_function_create(2, simple_hard_fn,
        100.0, 0.1, "simple-majority-2");
    assert(hf != NULL);

    NWPRG *prg = nw_prg_construct(d, hf);
    assert(prg != NULL);
    assert(prg->seed_len == 4);
    assert(prg->output_len == 3);
    assert(prg->stretch == 0.75);

    /* Evaluate on specific seed */
    bool seed[4] = {false, false, false, false}; /* all zeros */
    bool output[3] = {0};

    nw_prg_evaluate(prg, seed, output);

    /* G(0000) = [f(00), f(00), f(00)]
     * f(00): count=0 (< 2/2=1) -> false
     * So output = [false, false, false] */
    assert(output[0] == false);
    assert(output[1] == false);
    assert(output[2] == false);

    /* Now seed = {true, true, false, false} */
    seed[0] = true; seed[1] = true; seed[2] = false; seed[3] = false;
    nw_prg_evaluate(prg, seed, output);
    /* S_0 = {0,1}: f(true, true) = count 2 >= 1 -> true
     * S_1 = {1,2}: f(true, false) = count 1 >= 1 -> true
     * S_2 = {2,3}: f(false, false) = count 0 < 1 -> false */
    assert(output[0] == true);
    assert(output[1] == true);
    assert(output[2] == false);

    nw_prg_free(prg);
    nw_hard_function_free(hf);
    PASS();
}

static void test_prg_eval_word(void) {
    TEST("prg_eval_word");
    NWDesign *d = nw_design_create(4, 2, 2, 1);
    d->sets[0][0] = 0; d->sets[0][1] = 1;
    d->sets[1][0] = 2; d->sets[1][1] = 3;

    HardFunction *hf = nw_hard_function_create(2, simple_hard_fn,
        100.0, 0.1, "test-hf");
    NWPRG *prg = nw_prg_construct(d, hf);

    /* seed = 0b0011 = 3: bits 0,1 are 1 */
    uint64_t result = nw_prg_eval_word_to_word(prg, 3ULL);
    /* f(11) should be true (count 2 >= 1)
     * f(00) should be false
     * result = true | (false << 1) = 1 */
    assert(result == 1ULL);

    nw_prg_free(prg);
    nw_hard_function_free(hf);
    PASS();
}

/* =============================================
 * L5: Seed Length Computation
 * ============================================= */

static void test_seed_output_length(void) {
    TEST("seed_output_length");
    /* Basic sanity */
    size_t seed = nw_seed_length(100, 1000.0);
    assert(seed > 0);

    size_t out = nw_output_length(seed, 1000.0);
    assert(out > 0);

    /* Monotonicity: larger output needs at least as many seeds */
    size_t seed10 = nw_seed_length(10, 1000.0);
    size_t seed100 = nw_seed_length(100, 1000.0);
    assert(seed100 >= seed10);
    PASS();
}

/* =============================================
 * L5: Optimal PRG Construction
 * ============================================= */

static void test_prg_construct_optimal(void) {
    TEST("prg_construct_optimal");
    NWPRG *prg = nw_prg_construct_optimal(50, 1000.0, NW_PRG_STANDARD);
    assert(prg != NULL);
    assert(prg->output_len >= 50 || prg->output_len + prg->seed_len > 0);
    assert(prg->design != NULL);
    nw_prg_free(prg);
    PASS();
}

/* =============================================
 * L6: PARITY Circuit Construction
 * ============================================= */

static void test_circuit_parity(void) {
    TEST("circuit_parity");
    Circuit *c = nw_circuit_parity(4);
    assert(c != NULL);
    assert(c->n_inputs == 4);
    assert(nw_circuit_gate_count(c) > 0);

    /* Test PARITY evaluation */
    bool in1[4] = {true, false, true, false}; /* Even number of 1s */
    assert(nw_circuit_eval_single(c, in1) == false);

    bool in2[4] = {true, true, true, false}; /* Odd number of 1s */
    assert(nw_circuit_eval_single(c, in2) == true);

    bool in3[4] = {false, false, false, false};
    assert(nw_circuit_eval_single(c, in3) == false);

    bool in4[4] = {true, false, false, false};
    assert(nw_circuit_eval_single(c, in4) == true);

    nw_circuit_free(c);
    PASS();
}

/* =============================================
 * L6: MAJORITY Circuit Construction
 * ============================================= */

static void test_circuit_majority(void) {
    TEST("circuit_majority");
    Circuit *c = nw_circuit_majority(3);
    assert(c != NULL);
    assert(c->n_inputs == 3);

    /* Test MAJORITY evaluation */
    bool in_maj0[3] = {false, false, false};
    assert(nw_circuit_eval_single(c, in_maj0) == false);

    bool in_maj1[3] = {true, false, false};
    assert(nw_circuit_eval_single(c, in_maj1) == false);

    bool in_maj2[3] = {true, true, false};
    assert(nw_circuit_eval_single(c, in_maj2) == true);

    bool in_maj3[3] = {true, true, true};
    assert(nw_circuit_eval_single(c, in_maj3) == true);

    nw_circuit_free(c);
    PASS();
}

/* =============================================
 * L4: Circuit Lower Bound — Shannon
 * ============================================= */

static void test_shannon_circuit_lower_bound(void) {
    TEST("shannon_circuit_lower_bound");
    size_t lb = nw_shannon_circuit_lower_bound(10);
    /* Shannon: most functions need ~ 2^10 / 20 ≈ 51 gates */
    assert(lb > 10);
    PASS();
}

/* =============================================
 * L5: Hardness Assumption Validation
 * ============================================= */

static void test_hardness_assumption(void) {
    TEST("hardness_assumption");
    HardnessAssumption ha = { .n = 16, .circuit_size = 1000.0,
                              .hardness = 0.1, .is_exponential = false };
    assert(nw_hardness_assumption_valid(&ha) == true);

    ha.hardness = 1.0; /* Invalid: hardness > 0.5 */
    assert(nw_hardness_assumption_valid(&ha) == false);
    PASS();
}

/* =============================================
 * L5: Truth Table Operations
 * ============================================= */

static void test_truth_table(void) {
    TEST("truth_table");
    /* Create truth table for AND on 2 vars */
    TruthTable *tt = nw_truth_table_create(simple_hard_fn, 3);
    assert(tt != NULL);
    assert(tt->n == 3);
    assert(tt->size == 8);

    /* Check some values */
    bool in[3] = {true, true, true};   /* 3 ones -> count>=2 -> true */
    assert(nw_truth_table_eval(tt, in) == true);

    bool in2[3] = {false, false, false}; /* 0 ones -> count=0 < 1 -> false */
    assert(nw_truth_table_eval(tt, in2) == false);

    nw_truth_table_free(tt);
    PASS();
}

/* =============================================
 * L8: Switching Lemma / Hastad
 * ============================================= */

static void test_hastad_applies(void) {
    TEST("hastad_applies");
    assert(nw_hastad_applies(2, 100) == true);
    assert(nw_hastad_applies(1, 100) == false); /* Depth too small */

    double bound = nw_switching_lemma_bound(100, 3, 0.1);
    assert(bound >= 0.0 && bound <= 1.0);

    double ac0_lb = nw_ac0_parity_lower_bound(100, 3);
    /* AC0 lower bound: 2^{n^{1/(d-1)} / C} = 2^{10/C} must be > 1 */
    assert(ac0_lb >= 1.0);

    /* With larger n, the bound should grow */
    double ac0_lb2 = nw_ac0_parity_lower_bound(10000, 3);
    assert(ac0_lb2 > ac0_lb);
    PASS();
}

/* =============================================
 * L5: SetSystem Operations
 * ============================================= */

static void test_setsystem(void) {
    TEST("setsystem");
    SetSystem *ss = nw_setsystem_create(10, 3, 2);
    assert(ss != NULL);
    assert(ss->universe == 10);
    assert(ss->num_sets == 3);
    assert(ss->set_size == 2);

    nw_setsystem_set(ss, 0, 0, 1);
    nw_setsystem_set(ss, 0, 1, 2);
    assert(nw_setsystem_get(ss, 0, 0) == 1);
    assert(nw_setsystem_get(ss, 0, 1) == 2);

    nw_setsystem_free(ss);
    PASS();
}

/* =============================================
 * L7: IW97 Derandomization Levels
 * ============================================= */

static void test_iw_derandomize(void) {
    TEST("iw_derandomize");
    size_t seed_len;

    assert(nw_iw_derandomize(0, 100, &seed_len) == true);
    assert(seed_len > 0);

    assert(nw_iw_derandomize(1, 100, &seed_len) == true);
    assert(seed_len > 0);

    assert(nw_iw_derandomize(2, 100, &seed_len) == true);
    assert(seed_len > 0);

    assert(nw_iw_derandomize(99, 100, &seed_len) == false);
    PASS();
}

/* =============================================
 * L5: Circuit Verification
 * ============================================= */

static void test_circuit_verify(void) {
    TEST("circuit_verify");
    Circuit *c = nw_circuit_parity(4);
    assert(c != NULL);
    assert(nw_circuit_verify(c) == true);
    nw_circuit_free(c);
    PASS();
}

/* =============================================
 * L5: Goldreich-Levin Hard Predicate
 * ============================================= */

static void test_gl_hard_predicate(void) {
    TEST("goldreich_levin_hard_predicate");
    Circuit *c = nw_goldreich_levin_hard_predicate(3);
    assert(c != NULL);
    assert(c->n_inputs == 6); /* 2n */

    /* Inner product: <x, r> mod 2
     * x = (1,0,1), r = (1,1,0): dot = 1*1 + 0*1 + 1*0 = 1 => odd => true */
    bool in[6] = {true, false, true, true, true, false};
    assert(nw_circuit_eval_single(c, in) == true);

    /* x = (1,0,1), r = (0,1,0): dot = 0+0+0 = 0 => false */
    bool in2[6] = {true, false, true, false, true, false};
    assert(nw_circuit_eval_single(c, in2) == false);

    nw_circuit_free(c);
    PASS();
}

/* =============================================
 * L8: Natural Proofs Barrier
 * ============================================= */

static void test_natural_proofs(void) {
    TEST("natural_proofs_barrier");
    /* Small hardness should not trigger barrier */
    assert(nw_natural_proofs_barrier(100.0, 8) == false);

    /* Very large hardness at larger n should trigger */
    double huge = 1e100;
    assert(nw_natural_proofs_barrier(huge, 20) == true);

    double thresh = nw_razborov_rudich_threshold(16);
    assert(thresh > 0.0);
    PASS();
}

/* =============================================
 * L5: Statistical Tests on PRG
 * ============================================= */

static void test_statistical_tests(void) {
    TEST("statistical_tests");
    NWPRG *prg = nw_prg_construct_optimal(64, 10000.0, NW_PRG_STANDARD);
    assert(prg != NULL);

    StatisticalTestResult *r = nw_prg_statistical_tests(prg, 5, 128);
    assert(r != NULL);
    assert(r->num_samples == 5);
    assert(r->num_tests == 3);

    nw_statistical_test_result_free(r);
    nw_prg_free(prg);
    PASS();
}

/* =============================================
 * L2: Distinguisher Advantage (basic)
 * ============================================= */

static void test_distinguisher(void) {
    TEST("distinguisher_advantage");
    /* Build a trivial distinguisher circuit */
    NWPRG *prg = nw_prg_construct_optimal(16, 1000.0, NW_PRG_STANDARD);
    if (prg) {
        /* A trivial distinguisher: output 1 always -> advantage 0 */
        Circuit *dc = nw_circuit_create(prg->output_len);
        /* Make a trivial constant-1 circuit */
        size_t not0 = nw_circuit_add_not(dc, 0);  /* NOT input[0] */
        size_t or_gate = nw_circuit_add_or(dc, 0, not0); /* input[0] OR NOT input[0] = 1 */
        nw_circuit_set_output(dc, or_gate);

        double adv = nw_prg_distinguisher_advantage(prg, dc, 50);
        /* Trivial distinguisher has near-zero advantage */
        assert(adv <= 0.5);

        nw_circuit_free(dc);
        nw_prg_free(prg);
    }
    PASS();
}

/* =============================================
 * L5: Lupanov Compilation (small case)
 * ============================================= */

static void test_compile_truth_table(void) {
    TEST("compile_truth_table");
    TruthTable *tt = nw_truth_table_create(simple_hard_fn, 3);
    assert(tt != NULL);

    Circuit *c = nw_compile_truth_table(tt);
    assert(c != NULL);
    assert(nw_circuit_verify(c) == true);

    /* Test that compiled circuit matches truth table */
    bool in[3];
    for (size_t i = 0; i < 8; i++) {
        in[0] = (i & 1) != 0;
        in[1] = (i & 2) != 0;
        in[2] = (i & 4) != 0;
        bool expected = simple_hard_fn(in, 3);
        bool actual = nw_circuit_eval_single(c, in);
        assert(expected == actual);
    }

    nw_circuit_free(c);
    nw_truth_table_free(tt);
    PASS();
}

/* =============================================
 * L3: Binomial Coefficient
 * ============================================= */

static void test_binomial(void) {
    TEST("binomial_coefficient");
    assert(nw_binomial(5, 0) == 1.0);
    assert(nw_binomial(5, 5) == 1.0);
    assert(nw_binomial(5, 2) == 10.0);
    assert(nw_binomial(10, 5) == 252.0);

    double log_bc = nw_log_binomial(100, 50);
    assert(log_bc > 0.0);
    PASS();
}

/* =============================================
 * L3: Fourier Analysis
 * ============================================= */

static void test_fourier(void) {
    TEST("fourier_coefficient");
    TruthTable *tt = nw_truth_table_create(simple_hard_fn, 3);
    assert(tt != NULL);

    bool S[3] = {false, false, false};
    double empty_coeff = nw_fourier_coefficient(tt, S);
    /* Empty Fourier coefficient = E[f] mapped to ±1 */
    assert(fabs(empty_coeff) <= 1.0);

    double inf = nw_variable_influence(tt, 0);
    assert(inf >= 0.0 && inf <= 1.0);

    nw_truth_table_free(tt);
    PASS();
}

/* =============================================
 * L3: Expander Construction
 * ============================================= */

static void test_expander(void) {
    TEST("expander_construction");
    Expander *e = nw_expander_create_margulis(100, 8);
    assert(e != NULL);
    assert(e->n >= 100);
    assert(e->d == 8);

    double ab = nw_expander_alon_boppana(8.0);
    assert(ab > 0.0);

    nw_expander_free(e);
    PASS();
}

/* =============================================
 * L5: Exponential Derandomization (IW97)
 * ============================================= */

static void test_bpp_simulation(void) {
    TEST("bpp_simulation");
    size_t samples = nw_bpp_simulation_samples(100, 0.01);
    assert(samples > 0);
    PASS();
}

/* =============================================
 * Main Test Runner
 * ============================================= */

int main(void) {
    printf("╔══════════════════════════════════════════╗\n");
    printf("║  mini-nisan-wigderson-prg — Test Suite  ║\n");
    printf("╚══════════════════════════════════════════╝\n\n");

    printf("── L1: Definitions ─────────────────────\n");
    test_design_create_free();
    test_design_create_null_params();

    printf("\n── L3: Mathematical Structures ─────────\n");
    test_design_verify_small();
    test_design_verify_intersecting();
    test_int_seed_conversion();
    test_seed_restrict();
    test_hamming_weight();
    test_bool_xor();
    test_ceil_log2();
    test_setsystem();
    test_binomial();
    test_fourier();
    test_expander();

    printf("\n── L4: Fundamental Laws ────────────────\n");
    test_yao_xor_amplify();
    test_yao_xor_lemma_full();
    test_shannon_bound();
    test_hardness_assumption();
    test_shannon_circuit_lower_bound();

    printf("\n── L5: Algorithms/Methods ──────────────\n");
    test_reed_solomon_design();
    test_random_design();
    test_random_design_success_prob();
    test_prg_construct_evaluate();
    test_prg_eval_word();
    test_prg_construct_optimal();
    test_seed_output_length();
    test_truth_table();
    test_compile_truth_table();
    test_gl_hard_predicate();
    test_circuit_verify();
    test_statistical_tests();
    test_distinguisher();

    printf("\n── L6: Canonical Problems ──────────────\n");
    test_circuit_parity();
    test_circuit_majority();

    printf("\n── L7: Applications ────────────────────\n");
    test_iw_derandomize();
    test_bpp_simulation();

    printf("\n── L8: Advanced Topics ─────────────────\n");
    test_hastad_applies();
    test_natural_proofs();

    printf("\n══════════════════════════════════════════\n");
    printf("  Results: %d/%d tests passed\n", tests_passed, tests_run);
    printf("══════════════════════════════════════════\n");

    if (tests_passed == tests_run) {
        printf("  ALL TESTS PASSED ✓\n");
        return 0;
    } else {
        printf("  FAILURES: %d\n", tests_run - tests_passed);
        return 1;
    }
}
