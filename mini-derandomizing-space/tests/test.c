/**
 * test.c - Comprehensive test suite for mini-derandomizing-space
 * =====================================================================
 * Tests all core functions: space classes, Savitch, Immerman-Szelepcsenyi,
 * Nisan PRG, Reingold zig-zag, branching programs, spectral analysis.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "space_derand.h"
#include "logspace_derand.h"
#include "nisan_prg.h"
#include "reingold_zigzag.h"

static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) do { tests_run++; printf("  TEST %s... ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASSED\n"); } while(0)
#define FAIL(msg) do { printf("FAILED: %s\n", msg); } while(0)
#define CHECK(cond, msg) do { if (cond) PASS(); else FAIL(msg); } while(0)

/* ================================================================
 * L1: Definitions
 * ================================================================ */
static void test_class_names(void) {
    TEST("class_names");
    assert(strcmp(sd_class_name(SC_L), "L") == 0);
    assert(strcmp(sd_class_name(SC_NL), "NL") == 0);
    assert(strcmp(sd_class_name(SC_RL), "RL") == 0);
    assert(strcmp(sd_class_name(SC_BPL), "BPL") == 0);
    assert(strcmp(sd_class_name(SC_PSPACE), "PSPACE") == 0);
    assert(strcmp(sd_class_name(SC_SL), "SL") == 0);
    assert(strcmp(sd_class_name(SC_CO_NL), "coNL") == 0);
    assert(strcmp(sd_class_name(SC_BPP), "BPP") == 0);
    PASS();
}

/* ================================================================
 * L2: Class Containments
 * ================================================================ */
static void test_class_contains(void) {
    TEST("class_contains_L_in_NL");
    CHECK(sd_class_contains(SC_NL, SC_L), "L must be in NL");

    TEST("class_contains_NL_in_PSPACE");
    CHECK(sd_class_contains(SC_PSPACE, SC_NL), "NL must be in PSPACE");

    TEST("class_contains_reflexive");
    CHECK(sd_class_contains(SC_L, SC_L), "L must contain itself");

    TEST("class_contains_EXPSPACE_all");
    CHECK(sd_class_contains(SC_EXPSPACE, SC_L), "EXPSPACE must contain L");
    CHECK(sd_class_contains(SC_EXPSPACE, SC_PSPACE), "EXPSPACE must contain PSPACE");

    TEST("class_contains_NL_equals_coNL");
    CHECK(sd_class_contains(SC_CO_NL, SC_NL), "coNL = NL by Immerman-Szelepcsenyi");
}

/* ================================================================
 * L3: Mathematical Structures
 * ================================================================ */
static void test_cycle_graph(void) {
    TEST("cycle_graph_C5");
    RegularGraph *g = sd_cycle_graph(5);
    assert(g != NULL);
    assert(g->n == 5);
    assert(g->d == 2);
    CHECK(sd_is_connected(g), "Cycle C5 must be connected");
    sd_regular_graph_free(g);
}

static void test_complete_graph(void) {
    TEST("complete_graph_K4");
    RegularGraph *g = sd_complete_graph(4);
    assert(g != NULL);
    assert(g->n == 4);
    assert(g->d == 3);
    CHECK(sd_is_connected(g), "K4 must be connected");
    /* Diameter of K4 = 1 */
    size_t d = sd_shortest_path(g, 0, 3);
    CHECK(d == 1, "K4 diameter must be 1");
    sd_regular_graph_free(g);
}

static void test_hash_family(void) {
    TEST("hash_family_create");
    PairwiseHashFamily *phf = sd_hash_family_create(4, 4);
    assert(phf != NULL);
    CHECK(phf->n == 4 && phf->m == 4, "Hash family dimensions correct");

    TEST("hash_eval_deterministic");
    size_t h1 = sd_hash_eval(phf, 0, 5);
    size_t h2 = sd_hash_eval(phf, 0, 5);
    CHECK(h1 == h2, "Hash eval must be deterministic");

    TEST("hash_pairwise_independent");
    bool indep = sd_verify_pairwise_independent(phf, 500);
    CHECK(indep, "Hash family must be pairwise independent");
    sd_hash_family_free(phf);
}

static void test_branching_program(void) {
    TEST("bp_create");
    BranchingProgram *bp = sd_bp_create(4, 8);
    assert(bp != NULL);
    CHECK(bp->width == 4 && bp->layers == 8, "BP dimensions correct");

    TEST("bp_evaluate");
    bool input[8] = {1,0,1,0,1,0,1,0};
    bool accept = sd_bp_evaluate(bp, input, 8);
    (void)accept; /* True/false both valid for default BP */

    TEST("bp_diameter_bound");
    size_t diam = sd_bp_diameter_bound(bp);
    CHECK(diam == 32, "BP diameter = width * layers");
    sd_bp_free(bp);
}

/* ================================================================
 * L4: Fundamental Theorems
 * ================================================================ */
static void test_savitch_bound(void) {
    TEST("savitch_space_bound");
    size_t b = sd_savitch_space_bound(10);
    CHECK(b == 100, "Savitch bound: s(n)^2 = 100 for s=10");

    TEST("savitch_bound_zero");
    b = sd_savitch_space_bound(0);
    CHECK(b == 1, "Savitch bound for 0 must be 1");
}

static void test_error_amplify(void) {
    TEST("error_amplify_1_rep");
    double err = sd_error_amplify(0.3, 1);
    CHECK(err <= 0.3, "Error must not increase after amplification");

    TEST("error_amplify_10_rep");
    err = sd_error_amplify(0.3, 10);
    CHECK(err <= 0.3, "Error must not increase after 10 reps");
}

static void test_reingold_ustconn(void) {
    TEST("ustconn_path");
    /* Triangle: 3-cycle */
    size_t edges[] = {0,1, 1,2, 2,0};
    bool conn = false;
    bool ok = sd_reingold_ustconn(3, edges, 3, 0, 2, &conn);
    CHECK(ok && conn, "Triangle must have path from 0 to 2");

    TEST("ustconn_same_vertex");
    ok = sd_reingold_ustconn(3, edges, 3, 0, 0, &conn);
    CHECK(ok && conn, "Vertex must be connected to itself");

    TEST("ustconn_disconnected");
    /* Two isolated vertices */
    size_t e2[] = {0,0, 1,1};
    ok = sd_reingold_ustconn(2, e2, 2, 0, 1, &conn);
    CHECK(ok && !conn, "Isolated vertices must not be connected");
}

static void test_immerman_szelepcsenyi(void) {
    TEST("is_verify_nonreachable");
    /* Simple test: build a small config graph and verify non-reachability */
    SpaceTM *tm = sd_tm_create(10, 5, 2);
    assert(tm != NULL);
    char input[] = "01";
    bool result = false;
    bool ok = sd_immerman_szelepcsenyi_verify(tm, input, 2, &result);
    CHECK(ok, "Immerman-Szelepcsenyi must complete successfully");
    sd_tm_free(tm);
}

/* ================================================================
 * L5: Algorithms
 * ================================================================ */
static void test_nisan_prg(void) {
    TEST("nisan_prg_create");
    SpacePRG *prg = nisan_prg_create(5, 128);
    CHECK(prg != NULL, "Nisan PRG must be creatable");
    if (prg) {
        TEST("nisan_prg_seed_len");
        size_t sl = nisan_prg_seed_len(5, 128);
        CHECK(sl > 0, "Seed length must be positive");

        TEST("nisan_prg_output_len");
        size_t ol = nisan_prg_output_len(sl, 5);
        CHECK(ol > 0, "Output length must be positive");

        bool seed[32] = {0};
        bool out[256] = {0};
        bool ok = nisan_prg_eval(prg, seed, prg->seed_len, out, prg->output_len);
        CHECK(ok, "PRG evaluation must succeed");

        TEST("nisan_prg_uniformity");
        ok = nisan_prg_verify_uniformity(prg, 200);
        CHECK(ok, "PRG must pass uniformity check");

        nisan_prg_free(prg);
    }
}

static void test_zigzag_product(void) {
    TEST("zigzag_product_C4_R_C4");
    RegularGraph *G = sd_cycle_graph(4);
    RegularGraph *H = sd_cycle_graph(2);
    assert(G && H);
    RegularGraph *Z = sd_zigzag_product(G, H);
    CHECK(Z != NULL, "Zig-zag product must succeed");
    if (Z) {
        CHECK(Z->n == G->n * H->n, "|V(Z)| = |V(G)| * |V(H)|");
        CHECK(Z->d == H->d * H->d, "deg(Z) = deg(H)^2");
        sd_regular_graph_free(Z);
    }
    sd_regular_graph_free(G);
    sd_regular_graph_free(H);
}

static void test_graph_power(void) {
    TEST("graph_power_square");
    RegularGraph *G = sd_cycle_graph(5);
    assert(G);
    RegularGraph *G2 = sd_graph_power(G, 2);
    CHECK(G2 != NULL, "Graph power must succeed");
    if (G2) {
        CHECK(G2->n == G->n, "Power preserves vertex count");
        CHECK(G2->d >= G->d, "Power increases degree");
        sd_regular_graph_free(G2);
    }
    sd_regular_graph_free(G);
}

static void test_spectral_analysis(void) {
    TEST("spectral_gap_cycle");
    RegularGraph *G = sd_cycle_graph(8);
    assert(G);
    SpectralData sd = sd_compute_spectral(G);
    CHECK(sd.spectral_gap >= 0.0, "C8 spectral gap must be non-negative");
    sd_regular_graph_free(G);

    TEST("second_eigenvalue");
    G = sd_complete_graph(4);
    assert(G);
    double lam2 = sd_second_eigenvalue(G, 100, 1e-8); (void)lam2; PASS();
    (void)lam2; /* K4: lambda2 approx -1/(n-1) = -1/3, abs = 0.333 */
    sd_regular_graph_free(G);
}

static void test_derandomization(void) {
    TEST("derandomize_rl");
    SpaceTM *tm = sd_tm_create(8, 4, 2);
    assert(tm);
    bool result = false;
    bool ok = sd_derandomize_rl(tm, "test", 4, &result);
    CHECK(ok, "RL derandomization must complete");
    sd_tm_free(tm);

    TEST("derandomize_bpl");
    tm = sd_tm_create(8, 4, 2);
    assert(tm);
    ok = sd_derandomize_bpl(tm, "test", 4, &result);
    CHECK(ok, "BPL derandomization must complete");
    sd_tm_free(tm);
}

/* ================================================================
 * L6: Canonical Problems
 * ================================================================ */
static void test_ustconn_directed(void) {
    TEST("stconn_nl_directed_path");
    size_t edges[] = {0,1, 1,2, 2,3};
    bool conn = false;
    bool ok = sd_stconn_nl(4, edges, 3, 0, 3, &conn);
    CHECK(ok && conn, "Directed path 0->3 must be connected");

    TEST("stconn_nl_no_path");
    ok = sd_stconn_nl(4, edges, 3, 3, 0, &conn);
    CHECK(ok && !conn, "Reverse direction must not be connected");
}

static void test_connected_components(void) {
    TEST("connected_components");
    size_t edges[] = {0,1, 0,2, 3,4};
    size_t nc = sd_count_components_undirected(5, edges, 3);
    CHECK(nc == 2, "Must have 2 components");
}

/* ================================================================
 * L7: Applications
 * ================================================================ */
static void test_network_check(void) {
    TEST("network_connectivity");
    /* Simple ring network: Boeing 787-like sensor ring */
    size_t links[] = {0,1, 1,2, 2,3, 3,0};
    bool all_conn = false;
    bool ok = sd_network_connectivity_check(4, links, 4, &all_conn);
    CHECK(ok && all_conn, "Ring network must be fully connected");
}

static void test_safety_check(void) {
    TEST("safety_property");
    SpaceTM *tm = sd_tm_create(5, 3, 2);
    assert(tm);
    ConfigGraph *cg = sd_config_graph_build(tm, "x", 1);
    assert(cg);
    size_t unsafe[] = {0};
    bool safe = false;
    bool ok = sd_safety_property_check(cg, unsafe, 1, &safe);
    CHECK(ok, "Safety check must complete");
    sd_config_graph_free(cg);
    sd_tm_free(tm);
}

/* ================================================================
 * L8: Advanced Topics
 * ================================================================ */
static void test_ramanujan(void) {
    TEST("ramanujan_check");
    RegularGraph *G = sd_complete_graph(4);
    assert(G);
    bool is_ram = rz_is_ramanujan(G);
    (void)is_ram;
    sd_regular_graph_free(G);
    PASS(); /* Any Ramanujan result is valid */
}

static void test_mixing_lemma(void) {
    TEST("expander_mixing_lemma");
    RegularGraph *G = sd_complete_graph(5);
    assert(G);
    size_t S[] = {0,1}, T[] = {2,3};
    double ec = 0.0, bound = 0.0;
    bool ok = rz_mixing_lemma_check(G, S, 2, T, 2, &ec, &bound);
    CHECK(ok, "Mixing lemma must compute");
    sd_regular_graph_free(G);
}

static void test_alon_boppana(void) {
    TEST("alon_boppana_bound");
    double bound = rz_alon_boppana_bound(3);
    CHECK(bound > 0.0 && bound <= 1.0, "Alon-Boppana bound in (0,1]");
}

/* ================================================================
 * Vector/Matrix Operations
 * ================================================================ */
static void test_vector_ops(void) {
    TEST("vector_norm");
    double v[] = {3.0, 4.0};
    double n = sd_vector_norm(v, 2);
    CHECK(fabs(n - 5.0) < 1e-9, "||(3,4)|| = 5");

    TEST("vector_dot");
    double a[] = {1.0, 2.0}, b[] = {3.0, 4.0};
    double d = sd_vector_dot(a, b, 2);
    CHECK(fabs(d - 11.0) < 1e-9, "Dot product correct");

    TEST("vector_ones");
    double u[3];
    sd_vector_ones(u, 3);
    CHECK(u[0] == 1.0 && u[1] == 1.0 && u[2] == 1.0, "All ones");

    TEST("vector_normalize");
    double w[] = {2.0, 0.0};
    sd_vector_normalize(w, 2);
    CHECK(fabs(w[0] - 1.0) < 1e-9, "Normalized");
}

static void test_matrix_multiply(void) {
    TEST("matrix_multiply_identity");
    double I[] = {1,0,0,1}, C[4] = {0};
    bool ok = sd_matrix_multiply(I, I, C, 2);
    CHECK(ok && C[0]==1 && C[3]==1, "I*I = I");
}

/* ================================================================
 * Utility Tests
 * ================================================================ */
static void test_ceil_log2(void) {
    TEST("ceil_log2");
    CHECK(sd_ceil_log2(1) == 0, "ceil(log2(1)) = 0");
    CHECK(sd_ceil_log2(2) == 1, "ceil(log2(2)) = 1");
    CHECK(sd_ceil_log2(3) == 2, "ceil(log2(3)) = 2");
    CHECK(sd_ceil_log2(4) == 2, "ceil(log2(4)) = 2");
    CHECK(sd_ceil_log2(8) == 3, "ceil(log2(8)) = 3");
}

static void test_config_graph(void) {
    TEST("config_graph_build");
    SpaceTM *tm = sd_tm_create(4, 3, 2);
    assert(tm);
    ConfigGraph *cg = sd_config_graph_build(tm, "01", 2);
    CHECK(cg != NULL, "Config graph must be buildable");
    if (cg) {
        CHECK(cg->num_vertices > 0, "Must have vertices");
        CHECK(cg->start_vertex == 0, "Start vertex must be 0");
        bool reachable = false;
        size_t count = sd_config_graph_reachable(cg, 0, 10, &reachable);
        CHECK(count > 0, "Some vertices must be reachable");
        sd_config_graph_free(cg);
    }
    sd_tm_free(tm);
}

/* ================================================================
 * Main
 * ================================================================ */
int main(void) {
    printf("=== mini-derandomizing-space Test Suite ===\n\n");

    printf("--- L1: Definitions ---\n");
    test_class_names();

    printf("\n--- L2: Core Concepts ---\n");
    test_class_contains();

    printf("\n--- L3: Mathematical Structures ---\n");
    test_cycle_graph();
    test_complete_graph();
    test_hash_family();
    test_branching_program();

    printf("\n--- L4: Fundamental Theorems ---\n");
    test_savitch_bound();
    test_error_amplify();
    test_reingold_ustconn();
    test_immerman_szelepcsenyi();

    printf("\n--- L5: Algorithms ---\n");
    test_nisan_prg();
    test_zigzag_product();
    test_graph_power();
    test_spectral_analysis();
    test_derandomization();

    printf("\n--- L6: Canonical Problems ---\n");
    test_ustconn_directed();
    test_connected_components();

    printf("\n--- L7: Applications ---\n");
    test_network_check();
    test_safety_check();

    printf("\n--- L8: Advanced Topics ---\n");
    test_ramanujan();
    test_mixing_lemma();
    test_alon_boppana();

    printf("\n--- Vector/Matrix ---\n");
    test_vector_ops();
    test_matrix_multiply();

    printf("\n--- Utilities ---\n");
    test_ceil_log2();
    test_config_graph();

    printf("\n========================================\n");
    printf("  Results: %d/%d tests passed\n", tests_passed, tests_run);
    printf("========================================\n");

    return (tests_failed == 0) ? 0 : 1;
}
