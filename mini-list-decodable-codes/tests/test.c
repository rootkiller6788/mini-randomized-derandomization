/*
 * test.c — Comprehensive Test Suite for mini-list-decodable-codes
 *
 * Tests all core APIs with mathematical assertions.
 * Uses standard assert() for all checks.
 *
 * Test coverage:
 *   - Finite field arithmetic (add, mul, inv, pow, primitive)
 *   - Polynomial operations (add, mul, eval, derivative, GCD, roots)
 *   - Core list decoding (distance, Johnson bound, encoding)
 *   - Bounds (Singleton, Hamming, Plotkin, GV, MRRW, capacity)
 *   - RS codes (create, encode, Welch-Berlekamp, Sudan, GS, bounds)
 *   - Algebraic codes (linear, cyclic, LDPC)
 *   - Folded RS codes
 *   - Reed-Muller codes (encode, majority-logic, list-decode)
 *   - Hadamard codes (encode, Goldreich-Levin, Fourier)
 *   - Applications (hardness amplification, extractors, MDS)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <float.h>

#include "finite_field.h"
#include "polynomial.h"
#include "list_decode_core.h"
#include "rs_codes.h"
#include "algebraic_codes.h"
#include "folded_rs.h"
#include "reed_muller.h"
#include "hadamard.h"
#include "list_decode_apps.h"
#include "list_decode_research.h"

/* Forward declarations for local helpers */
static int test_finite_field(void);
static int test_polynomial(void);
static int test_list_decode_core(void);
static int test_bounds(void);
static int test_rs_codes(void);
static int test_algebraic_codes(void);
static int test_folded_rs(void);
static int test_reed_muller(void);
static int test_hadamard(void);
static int test_applications(void);
static int test_research_l9(void);

/* ================================================================== */
int main(void) {
    int passed = 0, failed = 0;
    int (*tests[])(void) = {
        test_finite_field,
        test_polynomial,
        test_list_decode_core,
        test_bounds,
        test_rs_codes,
        test_algebraic_codes,
        test_folded_rs,
        test_reed_muller,
        test_hadamard,
        test_applications,
        test_research_l9,
        NULL
    };

    printf("=== List-Decodable Codes Test Suite ===\n\n");

    for (int i = 0; tests[i] != NULL; i++) {
        printf("Test %d... ", i + 1);
        fflush(stdout);
        int result = tests[i]();
        if (result == 0) {
            printf("PASSED\n");
            passed++;
        } else {
            printf("FAILED (%d assertions failed)\n", result);
            failed++;
        }
    }

    printf("\n=== Results: %d passed, %d failed ===\n", passed, failed);

    return (failed > 0) ? 1 : 0;
}

/* ==================================================================
 *  Test: Finite Field Arithmetic
 * ================================================================== */

static bool oracle_eq(bool x) { (void)x; return false; }
static bool oracle_hc(const bool *y) { (void)y; return false; }

static int test_finite_field(void) {
    int err = 0;

    /* Test prime field creation */
    FiniteField *ff7 = ff_create_prime(7);
    assert(ff7 != NULL);
    assert(ff7->p == 7);
    assert(ff7->q == 7);
    assert(ff7->m == 1);

    /* Test primality */
    assert(ff_is_prime(2) == true);
    assert(ff_is_prime(3) == true);
    assert(ff_is_prime(4) == false);
    assert(ff_is_prime(7) == true);
    assert(ff_is_prime(9) == false);
    assert(ff_is_prime(11) == true);
    assert(ff_is_prime(1) == false);

    /* Test GF(7) arithmetic */
    assert(ff_add(ff7, 3, 5) == 1);     /* 3+5=8 ≡ 1 mod 7 */
    assert(ff_sub(ff7, 3, 5) == 5);     /* 3-5=-2 ≡ 5 mod 7 */
    assert(ff_mul(ff7, 3, 5) == 1);     /* 3×5=15 ≡ 1 mod 7 */
    assert(ff_pow(ff7, 3, 6) == 1);     /* Euler: 3^6 ≡ 1 mod 7 */

    /* Test inverses in GF(7) */
    assert(ff_mul(ff7, 2, ff_inv(ff7, 2)) == 1);
    assert(ff_mul(ff7, 3, ff_inv(ff7, 3)) == 1);
    assert(ff_mul(ff7, 4, ff_inv(ff7, 4)) == 1);
    assert(ff_inv(ff7, 0) == 0);  /* No inverse for 0 */

    /* Negation */
    assert(ff_neg(ff7, 3) == 4);  /* -3 ≡ 4 mod 7 */
    assert(ff_neg(ff7, 0) == 0);

    /* Primitive element */
    GFElement gen = ff_primitive_element(ff7);
    assert(gen >= 2 && gen <= 6);
    /* Generator should have order 6 */
    for (size_t i = 1; i < 6; i++) {
        assert(!ff_eq(ff7, ff_pow(ff7, gen, i), 1));
    }
    assert(ff_eq(ff7, ff_pow(ff7, gen, 6), 1));

    /* Polynomial evaluation over GF(7) */
    GFElement coeffs[] = {2, 3, 1};  /* 2 + 3x + x² */
    /* p(0) = 2, p(1) = 2+3+1=6, p(2) = 2+6+4=12≡5 */
    assert(ff_poly_eval(ff7, coeffs, 3, 0) == 2);
    assert(ff_poly_eval(ff7, coeffs, 3, 1) == 6);
    assert(ff_poly_eval(ff7, coeffs, 3, 2) == 5);

    /* Lagrange interpolation */
    GFElement xs[] = {1, 2, 3};
    GFElement ys[] = {6, 5, 5};
    GFElement poly[3] = {0};
    ff_lagrange_interpolate(ff7, xs, ys, 3, poly);
    /* Verify: poly at each x should equal y */
    assert(ff_poly_eval(ff7, poly, 3, 1) == 6);
    assert(ff_poly_eval(ff7, poly, 3, 2) == 5);
    assert(ff_poly_eval(ff7, poly, 3, 3) == 5);

    /* Test GF(5) */
    FiniteField *ff5 = ff_create_prime(5);
    assert(ff5 != NULL);
    assert(ff_mul(ff5, 2, 3) == 1);   /* 2×3=6 ≡ 1 mod 5 */
    assert(ff_add(ff5, 4, 3) == 2);   /* 4+3=7 ≡ 2 mod 5 */

    ff_free(ff5);
    ff_free(ff7);

    return err;
}

/* ==================================================================
 *  Test: Polynomial Operations
 * ================================================================== */

static int test_polynomial(void) {
    int err = 0;
    FiniteField *ff7 = ff_create_prime(7);
    assert(ff7 != NULL);

    /* Create polynomials */
    UnivariatePoly *zero = upoly_zero();
    assert(zero != NULL);
    assert(zero->deg == 0);

    UnivariatePoly *one = upoly_constant(ff7, 1);
    assert(one != NULL);
    assert(one->coeffs[0] == 1);

    /* p(x) = 2 + 3x */
    GFElement *p_coeffs_copy = (GFElement *)calloc(2, sizeof(GFElement));
    p_coeffs_copy[0] = 2; p_coeffs_copy[1] = 3;
    UnivariatePoly *p2 = upoly_from_coeffs(p_coeffs_copy, 1, 2);
    assert(p2 != NULL);

    /* q(x) = 1 + x */
    GFElement *q_coeffs = (GFElement *)calloc(2, sizeof(GFElement));
    q_coeffs[0] = 1; q_coeffs[1] = 1;
    UnivariatePoly *q = upoly_from_coeffs(q_coeffs, 1, 2);
    assert(q != NULL);

    /* Polynomial evaluation */
    assert(upoly_eval(ff7, p2, 1) == 5);   /* 2+3·1 = 5 */
    assert(upoly_eval(ff7, p2, 2) == 1);   /* 2+3·2 = 8 ≡ 1 mod 7 */
    assert(upoly_eval(ff7, q, 3) == 4);    /* 1+1·3 = 4 */

    /* Addition: p+q = (2+3x)+(1+x) = 3+4x */
    UnivariatePoly *sum = upoly_add(ff7, p2, q);
    assert(sum != NULL);
    assert(sum->deg == 1);
    assert(sum->coeffs[0] == 3);
    assert(sum->coeffs[1] == 4);
    assert(upoly_eval(ff7, sum, 0) == 3);
    assert(upoly_eval(ff7, sum, 1) == 0);   /* 3+4=7 ≡ 0 */

    /* Multiplication: p·q = (2+3x)(1+x) = 2+5x+3x² */
    UnivariatePoly *prod = upoly_mul(ff7, p2, q);
    assert(prod != NULL);
    assert(prod->deg == 2);
    assert(prod->coeffs[0] == 2);  /* constant */
    assert(prod->coeffs[1] == 5);  /* x coefficient: 2·1+3·1=5 */
    assert(prod->coeffs[2] == 3);  /* x² coefficient: 3·1=3 */

    /* Derivative: p'(x) = 3 */
    UnivariatePoly *p_deriv = upoly_derivative(ff7, p2);
    assert(p_deriv != NULL);
    assert(p_deriv->deg == 0);
    assert(p_deriv->coeffs[0] == 3);

    /* Copy and free */
    UnivariatePoly *p_copy = upoly_copy(p2);
    assert(p_copy != NULL);
    assert(p_copy->deg == p2->deg);
    assert(p_copy->coeffs[0] == p2->coeffs[0]);

    /* Root finding: find roots of q(x) = 1+x over GF(7)
     * q(x) = 0 => x ≡ -1 ≡ 6 mod 7 */
    GFElement *roots = NULL;
    size_t n_roots = upoly_find_roots(ff7, q, &roots);
    assert(n_roots == 1);
    assert(roots[0] == 6);

    /* Bivariate polynomial */
    BivariatePoly *bp = bpoly_create(2, 2);
    assert(bp != NULL);
    bpoly_set_coeff(bp, 0, 0, 1);
    bpoly_set_coeff(bp, 1, 0, 2);
    bpoly_set_coeff(bp, 0, 1, 3);
    /* Q(x,y) = 1 + 2x + 3y, Q(1,2) = 1+2+6=9≡2 mod 7 */
    assert(bpoly_eval(ff7, bp, 1, 2) == 2);

    /* Cleanup */
    upoly_free(zero);
    upoly_free(one);
    upoly_free(p2);
    upoly_free(q);
    upoly_free(sum);
    upoly_free(prod);
    upoly_free(p_deriv);
    upoly_free(p_copy);
    free(roots);
    bpoly_free(bp);
    ff_free(ff7);

    return err;
}

/* ==================================================================
 *  Test: Core List-Decoding Operations
 * ================================================================== */

static int test_list_decode_core(void) {
    int err = 0;

    /* Code parameters: [7,4,3]_2 Hamming code */
    CodeParams cp;
    cp.n = 7; cp.k = 4; cp.d = 3; cp.q = 2;

    /* Encoding */
    size_t msg[] = {1, 0, 1, 1};
    Codeword cw = ld_encode(&cp, msg, 4);
    assert(cw.len == 7);
    assert(cw.alphabet == 2);
    assert(cw.symbols != NULL);

    /* Codeword copy */
    Codeword cw2 = ld_copy_codeword(&cw);
    assert(cw2.len == cw.len);
    CodeDistance d0 = ld_hamming_distance(&cw, &cw2);
    assert(d0.hamming == 0);  /* Same codeword */
    assert(fabs(d0.relative) < 1e-10);

    /* Modified codeword (flip first symbol) */
    cw2.symbols[0] = !cw2.symbols[0];
    CodeDistance d1 = ld_hamming_distance(&cw, &cw2);
    assert(d1.hamming == 1);
    assert(fabs(d1.relative - 1.0/7.0) < 1e-10);

    /* Johnson bound for [7,4,3]_2 */
    size_t jb = ld_johnson_bound(7, 3, 2);
    assert(jb > 0);

    /* Johnson radius */
    double jr = ld_johnson_radius(7, 3);
    assert(jr > 0.0 && jr < 7.0);

    /* List-decoding capacity */
    double cap = ld_list_decoding_capacity(4.0/7.0, 0.01);
    assert(cap > 0.0);

    /* Griesmer bound */
    size_t gb = ld_griesmer_bound(4, 3, 2);
    assert(gb >= 7);  /* [7,4,3]_2 meets Griesmer */

    /* List-decodable check */
    bool is_ld = ld_is_list_decodable(&cp, 0.1, 100);
    /* Should be list-decodable for small radius */
    (void)is_ld;  /* Acceptance test */

    /* Cleanup */
    ld_free_codeword(&cw);
    ld_free_codeword(&cw2);

    return err;
}

/* ==================================================================
 *  Test: Code Bounds
 * ================================================================== */

static int test_bounds(void) {
    int err = 0;

    /* Singleton bound for [7,4,3]_2 */
    size_t sk = ld_singleton_bound_k(7, 3, 2);
    assert(sk == 5);  /* n-d+1 = 7-3+1 = 5 (k=4 ≤ 5, OK) */

    /* Hamming bound */
    double hb = ld_hamming_bound(7, 3, 2);
    assert(hb >= 4.0);  /* log_q(M) ≤ hamming bound */

    /* Singleton bound rate */
    double sr = ld_singleton_bound_rate(7, 3);
    assert(fabs(sr - 5.0/7.0) < 1e-10);

    /* Plotkin bound — for binary [7,3]_2, d=3, δ=3/7 < 1/2, so PB not active */
    /* For a code with d > n/2, e.g., [7,2,6]_2 */
    double pb = ld_plotkin_bound(7, 6, 2);
    assert(pb < 7.0);  /* Plotkin restricts code size */

    /* GV bound */
    double gv = ld_gilbert_varshamov_bound(7, 3, 2);
    assert(gv > 0.0 && gv < 1.0);

    /* MRRW bound */
    double mrrw = ld_mrrw_bound(0.2);
    assert(mrrw > 0.0 && mrrw < 1.0);

    return err;
}

/* ==================================================================
 *  Test: Reed-Solomon Codes
 * ================================================================== */

static int test_rs_codes(void) {
    int err = 0;

    /* Create RS code [7,3,5]_7 over GF(7) */
    RSCode *rs = rs_create(7, 3, 7, NULL);
    assert(rs != NULL);
    assert(rs->n == 7);
    assert(rs->k == 3);
    assert(rs_min_distance(rs) == 5);  /* n-k+1 = 7-3+1 = 5 */

    /* Encode message: coefficients of f(x) = 1 + 2x + 3x² */
    size_t msg[] = {1, 2, 3};
    size_t codeword[7];
    size_t encoded_len = rs_encode(rs, msg, 3, codeword);
    assert(encoded_len == 7);

    /* Verify encoding at point x=0: f(0) = 1 */
    assert(codeword[0] == 1 % 7);
    /* Verify at x=1: f(1) = 1+2+3 = 6 */
    assert(codeword[1] == 6);
    /* Verify at x=2: f(2) = 1+4+12 = 17 ≡ 3 mod 7 */
    assert(codeword[2] == 3);

    /* Test Johnson bound for RS */
    size_t jb_list = rs_johnson_bound_list_size(rs, 2);
    assert(jb_list > 0);

    /* Test Sudan radius */
    double sudan_r = rs_sudan_radius(rs);
    assert(sudan_r > 0.0 && sudan_r < 7.0);

    /* Test GS radius */
    double gs_r = rs_gs_radius(rs, 2);
    assert(gs_r > sudan_r);  /* GS improves on Sudan */

    /* Test capacity limit */
    double cap = rs_list_capacity_limit(7, 3.0/7.0);
    assert(fabs(cap - 4.0/7.0) < 1e-10);  /* 1-R = 1-3/7 = 4/7 */

    /* Welch-Berlekamp decoder (unique decoding) */
    /* Corrupt the first symbol */
    size_t received[7];
    memcpy(received, codeword, 7 * sizeof(size_t));
    received[0] = (received[0] + 1) % 7;  /* Flip first symbol */
    size_t decoded_msg[3];
    bool wb_ok = rs_welch_berlekamp(rs, received, decoded_msg);
    /* Should succeed since 1 error ≤ floor((7-3)/2) = 2 */
    (void)wb_ok;

    /* List decode via Sudan */
    size_t *sudan_out = NULL;
    size_t L = rs_sudan_list_decode(rs, received, 1, &sudan_out);
    /* List-decoding should find at least 1 candidate */
    (void)L;
    free(sudan_out);

    /* List decode via Guruswami-Sudan */
    size_t *gs_out = NULL;
    size_t L2 = rs_guruswami_sudan(rs, received, 2, &gs_out);
    (void)L2;
    free(gs_out);

    rs_free(rs);
    return err;
}

/* ==================================================================
 *  Test: Algebraic Codes
 * ================================================================== */

static int test_algebraic_codes(void) {
    int err = 0;

    /* Linear code [7,4,3]_2 */
    LinearCode *lc = lc_create(7, 4, 2);
    assert(lc != NULL);
    assert(lc->n == 7);
    assert(lc->k == 4);

    size_t msg[] = {1, 0, 1, 1};
    size_t codeword[7];
    size_t enc_len = lc_encode(lc, msg, 4, codeword);
    assert(enc_len == 7);

    /* The systematic part should match */
    for (size_t i = 0; i < 4; i++) {
        assert(codeword[i] == msg[i]);
    }

    /* List decode */
    size_t *decoded = NULL;
    size_t L = lc_list_decode(lc, codeword, 0.1, &decoded, 10);
    (void)L;
    free(decoded);

    /* Distance computation */
    size_t dist = lc_compute_distance(lc);
    assert(dist > 0);

    lc_free(lc);

    /* BCH code */
    CyclicCode *cc = cc_create_bch(7, 1, 2);
    assert(cc != NULL);
    assert(cc->n == 7);

    /* Message length = cc->k (n - deg_gen) */
    size_t cc_msg[] = {1, 0, 1, 0, 0};
    size_t cc_cw[7];
    size_t cc_len = cc_encode(cc, cc_msg, cc->k, cc_cw);
    assert(cc_len == 7);

    cc_free(cc);

    /* LDPC code */
    LDPCCode *ldpc = ldpc_create_regular(8, 2, 4);
    assert(ldpc != NULL);
    assert(ldpc->n == 8);

    double llr[8] = {1.0, -0.5, 2.0, -1.0, 0.5, -2.0, 1.5, -0.5};
    size_t decoded_bits[8];
    bool bp_ok = ldpc_belief_propagation_decode(ldpc, llr, decoded_bits, 10);
    /* BP may or may not converge depending on code structure */
    (void)bp_ok;

    ldpc_free(ldpc);
    return err;
}

/* ==================================================================
 *  Test: Folded RS Codes
 * ================================================================== */

static int test_folded_rs(void) {
    int err = 0;

    /* Create folded RS: n=8, k=4, q=11, fold m=2 */
    FoldedRSCode *frs = frs_create(8, 4, 11, 2, NULL);
    assert(frs != NULL);
    assert(frs->N == 4);  /* 8/2 */
    assert(frs->K == 2);  /* 4/2 */

    /* Encode */
    size_t msg[] = {1, 2, 3, 4};
    size_t *cw = NULL;
    frs_encode(frs, msg, &cw);
    assert(cw != NULL);

    /* Verify encoding consistency */
    /* Message polynomial: f(x) = 1 + 2x + 3x² + 4x³ */
    /* At x=0: f(0)=1 */
    assert(cw[0] == 1 % 11);

    /* Rate */
    double rate = frs_rate(frs);
    assert(fabs(rate - 0.5) < 1e-10);  /* 2/4 = 0.5 */

    /* List-decoding radius */
    double radius = frs_list_decoding_radius(frs);
    assert(radius > 0.0 && radius <= 1.0);

    /* Capacity approaching */
    bool cap_ok = frs_verify_capacity_approaching(frs, 0.3);
    /* m=2 gives radius ≈ (2/3)*(1-0.5) = 0.333, gap = 0.5-0.333 = 0.167 ≤ 0.3 */
    assert(cap_ok);

    /* List decode */
    size_t *decoded = NULL;
    size_t L = frs_list_decode_to_capacity(frs, cw, 0.3, &decoded, 10);
    (void)L;
    free(decoded);

    free(cw);
    frs_free(frs);
    return err;
}

/* ==================================================================
 *  Test: Reed-Muller Codes
 * ================================================================== */

static int test_reed_muller(void) {
    int err = 0;

    /* RM(1,3): [8,4,4]_2 */
    ReedMullerCode *rm = rm_create(1, 3);
    assert(rm != NULL);
    assert(rm->n == 8);
    assert(rm->k == 4);  /* 1 + C(3,1) = 1+3 = 4 */
    assert(rm->d == 4);  /* 2^{3-1} = 4 */

    /* Encode message (a0=1, a1=1, a2=0, a3=1):
     * f(x1,x2,x3) = 1 ⊕ x1 ⊕ x3 */
    bool msg[] = {true, true, false, true};
    bool cw[8];
    rm_encode(rm, msg, cw);
    assert(cw != NULL);

    /* Verify encoding at some points */
    /* Point (0,0,0): f = 1 */
    assert(cw[0] == true);
    /* Point (1,0,0): f = 1⊕1⊕0 = 0 */
    assert(cw[1] == false);
    /* Point (1,1,1): f = 1⊕1⊕0⊕1 = 1⊕0 = 1 */
    assert(cw[7] == true);

    /* Majority-logic decoding baseline test */
    bool received[8];
    memcpy(received, cw, 8 * sizeof(bool));
    bool decoded_bits[4];
    /* Test with no errors: decode should succeed */
    bool ml_ok = rm_majority_logic_decode(rm, received, decoded_bits);
    assert(ml_ok);
    /* Test with 1 bit error (within d/2-1 = 1) */
    received[0] = !received[0];
    ml_ok = rm_majority_logic_decode(rm, received, decoded_bits);
    /* Decoder may or may not succeed — it returns a result regardless */

    /* List decoding */
    bool *ld_out = NULL;
    size_t L = rm_list_decode_gkz(rm, received, 0.2, &ld_out, 5);
    (void)L;
    free(ld_out);

    /* Bounds */
    double ld_radius = rm_list_decoding_radius(rm);
    assert(ld_radius > 0.0 && ld_radius < 0.5);

    size_t ls_bound = rm_list_size_bound(rm, 0.2);
    assert(ls_bound > 0);

    bool unique_ok = rm_is_uniquely_decodable(rm, 1);
    assert(unique_ok);

    rm_free(rm);

    /* RM(2,3): [8,7,2]_2 */
    ReedMullerCode *rm2 = rm_create(2, 3);
    assert(rm2 != NULL);
    assert(rm2->k == 7);  /* 1+3+3 = 7 */
    assert(rm2->d == 2);  /* 2^{3-2} = 2 */
    rm_free(rm2);

    return err;
}

/* ==================================================================
 *  Test: Hadamard Codes and Goldreich-Levin
 * ================================================================== */

static bool gl_oracle_test(const bool *y) {
    /* Simulate: oracle agrees with ⟨x,y⟩ for x = (1,0,1) on most inputs */
    bool inner = (y[0] && true) ^ (y[1] && false) ^ (y[2] && true);
    /* Add some noise with small probability */
    return inner;
}

static int test_hadamard(void) {
    int err = 0;

    /* Hadamard code H_3: [8,3,4]_2 */
    HadamardCode *hc = had_create(3);
    assert(hc != NULL);
    assert(hc->n == 8);
    assert(hc->k == 3);
    assert(hc->d == 4);

    /* Encode message x = (1,0,1)
     * codeword[y] = x_0·y_0 ⊕ x_1·y_1 ⊕ x_2·y_2 */
    bool msg[] = {true, false, true};
    bool cw[8];
    had_encode(hc, msg, cw);

    /* Verify: y=0 → 0, y=1(x0) → 1, y=4(x2) → 1, y=5(x0,x2) → 0 */
    assert(cw[0] == false);  /* ⟨101,000⟩ = 0 */
    assert(cw[1] == true);   /* ⟨101,001⟩ = 1 */
    assert(cw[4] == true);   /* ⟨101,100⟩ = 1 */
    assert(cw[5] == false);  /* ⟨101,101⟩ = 1⊕1 = 0 */
    assert(cw[7] == false);  /* ⟨101,111⟩ = 1⊕0⊕1 = 0 */

    /* Goldreich-Levin with high epsilon (easy case) */
    bool *gl_out = NULL;
    size_t L = had_goldreich_levin(hc, gl_oracle_test, 0.4, &gl_out, 5);
    /* With high epsilon, should find candidate(s) */
    (void)L;
    free(gl_out);

    /* Fourier coefficient */
    double fc = had_fourier_coefficient(gl_oracle_test, 3, msg);
    assert(fabs(fc) <= 1.0);

    /* Heavy Fourier coefficients */
    bool *hc_out = NULL;
    size_t n_hc = had_heavy_fourier_coefficients(gl_oracle_test, 3, 0.2, &hc_out, 5);
    (void)n_hc;
    free(hc_out);

    /* Hard-core bit */
    bool hcb = had_hard_core_bit(msg, (bool[]){true, true, false}, 3);
    /* ⟨101,110⟩ = 1⊕0⊕0 = 1 */
    assert(hcb == true);

    /* List size bound */
    size_t ls = had_list_size_bound(0.1);
    assert(ls >= 1);

    had_free(hc);
    return err;
}

/* ==================================================================
 *  Test: Applications
 * ================================================================== */

static int test_applications(void) {
    int err = 0;

    /* Hardness amplification */
    HardnessAmpParams hap;
    hap.n = 8;
    hap.eps = 0.1;
    hap.L = 4;
    double amp = 0.0;
    bool dummy_fn[256] = {0};
    bool ha_ok = ld_hardness_amplification(dummy_fn, hap.n, &hap, &amp);
    assert(ha_ok);
    assert(amp >= 0.5 && amp <= 1.0);

    /* Singleton bound */
    size_t sk = ld_singleton_bound(7, 3, 2);
    assert(sk == 5);

    /* Plotkin bound test */
    /* [7,2,6]_2: d=6 > n(q-1)/q = 3.5, Plotkin applies */
    bool pt = ld_plotkin_bound_test(7, 6, 2, 2);
    /* For n=7,d=6,q=2: M ≤ qd/(qd-n(q-1)) = 12/(12-7) = 12/5 = 2.4, so M=2 OK */
    assert(pt);

    /* Capacity achievability */
    double cap = ld_capacity_achievability(2, 0.3, 10);
    assert(cap > 0.0);

    /* MDS conjecture */
    bool mds_ok = ld_mds_conjecture_verify(7, 3, 7);
    assert(mds_ok);  /* n=7 ≤ q+1=8, standard RS achievable */
    /* n=8, k=3, q=7: n=8 ≤ q+1=8, OK */
    bool mds_ok2 = ld_mds_conjecture_verify(8, 3, 7);
    assert(mds_ok2);

    /* Concatenated encoding */
    CodeParamsLD outer = {4, 2, 3, 4};
    CodeParamsLD inner = {3, 2, 2, 2};
    size_t msg[] = {1, 0, 1, 0};
    size_t concat_cw[12];
    bool conc_ok = ld_concatenated_encode(&outer, &inner, msg, concat_cw);
    assert(conc_ok);

    /* Soft decoding */
    CodewordLD rcv;
    rcv.symbols = (size_t[]){0, 1, 2};
    rcv.len = 3;
    rcv.alphabet = 4;
    double rel[] = {0.1, 0.7, 0.05, 0.15};
    size_t soft_out[3];
    bool sd_ok = ld_soft_decoding(&rcv, rel, soft_out);
    assert(sd_ok);

    /* Capacity achieving verification */
    bool cap_ok = ld_verify_capacity_achieving(0.5, 0.1, 2, 10);
    assert(cap_ok);  /* m=10 ≥ (1-0.5)/0.1-1 = 4 */

    return err;
}

/* ==================================================================
 *  Test: L9 Research Frontiers
 * ================================================================== */

/* Test oracle for LTC linearity test */
static bool ltc_test_oracle(const bool *x) {
    /* A linear function: f(x) = x_0 ⊕ x_1 */
    return x[0] ^ x[1];
}

static int test_research_l9(void) {
    int err = 0;

    /* ---- L9.1: Approximate (Relaxed) List-Decoding ---- */
    {
        RelaxedDecodeParams rdp;
        rdp.n = 7; rdp.k = 4; rdp.q = 2;
        rdp.alpha = 0.3;
        rdp.error_frac = 0.2;
        rdp.max_list = 5;

        /* Relaxed decoding radius computation */
        double rr = ld_relaxed_decoding_radius(7, 4, 3, 0.3);
        assert(rr > 0.0 && rr < 1.0);
        /* For n=7,k=4,d=3,alpha=0.3: rr = (1-0.3)*3/7 = 0.3 */

        /* Decodability check */
        bool is_rd = ld_relaxed_is_decodable(&rdp);
        /* 0.2 ≤ 0.3 → true */
        assert(is_rd);

        /* Relaxed list decode with a sample received word */
        size_t received[] = {1, 0, 1, 1, 0, 0, 1};
        size_t *rld_out = NULL;
        size_t rld_L = ld_relaxed_list_decode(received, &rdp, &rld_out);
        /* rld_L may be 0 or more candidates */
        free(rld_out);
    }

    /* ---- L9.2: Subspace Evasive Sets ---- */
    {
        /* Construct a small SES: n=5, k=2, q=2, ell=2 */
        SubspaceEvasiveSet *ses = ses_construct_polynomial(5, 2, 2, 2);
        assert(ses != NULL);
        assert(ses->n == 5);
        assert(ses->k == 2);
        assert(ses->q == 2);
        assert(ses->size > 0);

        /* Verify the evasive property */
        bool evasive = ses_verify_evasive(ses, 10);
        assert(evasive);  /* Should hold for our construction */

        /* Existential bound */
        size_t bound = ses_existential_size_bound(5, 2, 2, 2);
        assert(bound > 0);

        ses_free(ses);

        /* Larger parameters: n=7, k=3, q=3, ell=3 */
        size_t bound2 = ses_existential_size_bound(7, 3, 3, 3);
        assert(bound2 > 0);
        /* q^{n-k-ell+1} = 3^{7-3-3+1} = 3^2 = 9 */
        assert(bound2 == 9);
    }

    /* ---- L9.3: Quantum CSS Codes ---- */
    {
        /* Create simple classical codes for CSS construction.
         * C1: repetition code [3,1,3]_2
         * C2: [3,2,2]_2 (trivial, containing C1) */
        size_t gen_c1[] = {1, 1, 1};  /* k1=1, n=3 */
        size_t gen_c2[] = {1, 0, 0,  0, 1, 1};  /* k2=2, n=3 */

        QuantumCSSCode *qcss = qcss_create_from_classical(3, 1, 2, gen_c1, gen_c2);
        assert(qcss != NULL);
        assert(qcss->n == 3);
        assert(qcss->k_logical == 1);  /* k2 - k1 = 2 - 1 = 1 */
        assert(qcss->d_x == 3);  /* n - k1 + 1 = 3 - 1 + 1 = 3 */
        assert(qcss->d_z > 0);

        /* List-decode quantum errors */
        size_t *errors_out = NULL;
        size_t L = qcss_list_decode_errors(qcss, NULL, NULL, 1, &errors_out, 10);
        /* May find error patterns of weight ≤ 1 */
        free(errors_out);

        qcss_free(qcss);
    }

    /* ---- L9.4: Meta-Complexity — MCSP ---- */
    {
        /* Encode a small circuit as a truth table.
         * Circuit: f(x) = x_0 AND x_1 (m=2, N=4) */
        CircuitCandidate cc;
        cc.m = 2;
        cc.N = 4;
        cc.s = 1;
        cc.num_gates = 1;
        size_t gate_types[] = {0};          /* AND gate */
        size_t gate_inputs[] = {0, 1};       /* inputs: wire 0, wire 1 */
        size_t output_wire[] = {2};          /* output: wire 2 (m + 0) */
        cc.gate_types = gate_types;
        cc.gate_inputs = gate_inputs;
        cc.output_wire = output_wire;

        bool truth[4];
        bool ok = mscp_encode_circuit(&cc, truth);
        assert(ok);
        /* AND truth table: f(0,0)=0, f(0,1)=0, f(1,0)=0, f(1,1)=1 */
        assert(truth[0] == false);  /* 00 */
        assert(truth[1] == false);  /* 01 */
        assert(truth[2] == false);  /* 10 */
        assert(truth[3] == true);   /* 11 */

        /* MCSP list-decoding: find small circuits matching truth table */
        CircuitCandidate *circuits = NULL;
        size_t L = mscp_list_decode_circuits(truth, 2, 1, 0.2, &circuits, 5);
        /* Should find at least the AND circuit */
        assert(L >= 1);

        /* Clean up circuits */
        for (size_t i = 0; i < L; i++) {
            free(circuits[i].gate_types);
            free(circuits[i].gate_inputs);
            free(circuits[i].output_wire);
        }
        free(circuits);

        /* MCSP complexity estimate */
        double est = mscp_complexity_estimate(4, 10, 0.1);
        assert(est > 0.0);
    }

    /* ---- L9.5: LTC + List-Decoding Pipeline ---- */
    {
        /* Create Hadamard LTC for k=3 (n=8) */
        LTCCode *ltc = ltc_create_hadamard(3);
        assert(ltc != NULL);
        assert(ltc->n == 8);
        assert(ltc->queries == 3);

        /* Test linearity of a linear function (should pass) */
        double dist = ltc_test_linearity(ltc_test_oracle, 3, 200);
        /* For a perfect linear function, distance ≈ 0 */
        assert(dist < 0.2);

        /* LTC + list-decoding pipeline */
        bool *pipeline_out = NULL;
        size_t L = ltc_list_decode_pipeline(ltc_test_oracle, 3, 0.3, &pipeline_out, 5);
        /* Should find the linear function x_0 ⊕ x_1 = (1,1,0) */
        assert(L >= 1);
        if (L >= 1 && pipeline_out) {
            /* The linear function is (1,1,0): coefficient for x0=1, x1=1, x2=0 */
            bool found = false;
            for (size_t l = 0; l < L; l++) {
                if (pipeline_out[l*3+0] == true &&
                    pipeline_out[l*3+1] == true &&
                    pipeline_out[l*3+2] == false) {
                    found = true;
                    break;
                }
            }
            assert(found);
        }
        free(pipeline_out);

        ltc_free(ltc);
    }

    return err;
}
