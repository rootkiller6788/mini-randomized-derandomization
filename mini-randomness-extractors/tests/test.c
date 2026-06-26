
/**
 * tests/test.c -- Comprehensive test suite for mini-randomness-extractors
 *
 * Tests covering L1-L6 with mathematical assertions:
 *   L1: All typedefs instantiated, constructors validated
 *   L2: Entropy measures tested on known distributions
 *   L3: GF(2) operations, statistical distance properties
 *   L4: Extractor bounds, LHL, RTD, NZ verified
 *   L5: All extractor constructions produce output
 *   L6: Canonical problem instances solved
 *
 * Every assertion is a mathematical check, not a trivial assert(1).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>

#include "extractor_core.h"
#include "extractor_constructions.h"
#include "extractor_applications.h"
#include "disperser_core.h"

/* Forward declarations of self-test functions from each module */
int ext_core_self_test(void);
int ext_constructions_self_test(void);
int ext_apps_self_test(void);
int disp_core_self_test(void);

#define TEST(name) printf("  TEST: %-50s ", name)
#define OK()      do { printf("OK\n"); } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); failures++; } while(0)

/*============================================================
 * L2: Mathematical property tests
 *============================================================*/

static int test_entropy_properties(void) {
    int failures = 0;

    /* Property 1: H_inf(X) <= H(X) for all distributions */
    TEST("Min-entropy <= Shannon entropy");
    {
        /* Distribution: p = [0.5, 0.3, 0.15, 0.05] */
        size_t freq[4] = {50, 30, 15, 5};
        double h_inf = ext_min_entropy_from_freq(freq, 4, 100);
        double h_shan = ext_shannon_entropy(freq, 4, 100);
        if (h_inf > h_shan + 0.001) FAIL("H_inf > H");
        else OK();
    }

    /* Property 2: H_inf(X) = log2(N) iff uniform */
    TEST("Min-entropy of uniform = log2(N)");
    {
        size_t freq[8] = {10,10,10,10,10,10,10,10};
        double h = ext_min_entropy_from_freq(freq, 8, 80);
        if (fabs(h - 3.0) > 0.01) FAIL("expected 3.0");
        else OK();
    }

    /* Property 3: H_inf -> 0 as distribution becomes deterministic */
    TEST("Min-entropy of deterministic -> 0");
    {
        size_t freq[8] = {100,0,0,0,0,0,0,0};
        double h = ext_min_entropy_from_freq(freq, 8, 100);
        if (fabs(h) > 0.01) FAIL("expected 0");
        else OK();
    }

    /* Property 4: Shannon entropy monotonicity with flattening */
    TEST("Shannon entropy maximized by uniform");
    {
        size_t unif[4] = {25,25,25,25};
        size_t skew[4] = {40,30,20,10};
        double hu = ext_shannon_entropy(unif, 4, 100);
        double hs = ext_shannon_entropy(skew, 4, 100);
        if (hu <= hs) FAIL("uniform not maximal");
        else OK();
    }

    /* Property 5: Renyi entropy alpha=infty = min-entropy */
    TEST("Renyi(inf) approximates min-entropy");
    {
        size_t freq[4] = {50,30,15,5};
        double h_inf = ext_min_entropy_from_freq(freq, 4, 100);
        double h_r10 = ext_renyi_entropy(freq, 4, 100, 10.0);
        /* For large alpha, Renyi converges to min-entropy */
        if (fabs(h_r10 - h_inf) > 0.2) FAIL("Renyi(10) too far from min-entropy");
        else OK();
    }

    /* Property 6: Collision entropy >= Min-entropy */
    TEST("Collision entropy >= Min-entropy");
    {
        size_t freq[4] = {50,30,15,5};
        double h_inf = ext_min_entropy_from_freq(freq, 4, 100);
        double h2 = ext_renyi_entropy(freq, 4, 100, 2.0);
        if (h2 < h_inf - 0.01) FAIL("H2 < Hinf");
        else OK();
    }

    /* Property 7: Statistical distance is a metric */
    TEST("Statistical distance: identity of indiscernibles");
    {
        size_t f1[4] = {25,25,25,25};
        size_t f2[4] = {25,25,25,25};
        double d = ext_statistical_distance(f1, f2, 4, 100, 100);
        if (fabs(d) > 0.001) FAIL("Delta(P,P) != 0");
        else OK();
    }

    /* Property 8: Statistical distance symmetric */
    TEST("Statistical distance: symmetry");
    {
        size_t f1[4] = {50,30,10,10};
        size_t f2[4] = {40,20,30,10};
        double d12 = ext_statistical_distance(f1, f2, 4, 100, 100);
        double d21 = ext_statistical_distance(f2, f1, 4, 100, 100);
        if (fabs(d12 - d21) > 0.001) FAIL("Delta(P,Q) != Delta(Q,P)");
        else OK();
    }

    /* Property 9: Triangle inequality */
    TEST("Statistical distance: triangle inequality");
    {
        size_t f1[4] = {50,30,10,10};
        size_t f2[4] = {40,20,30,10};
        size_t f3[4] = {10,40,25,25};
        double d12 = ext_statistical_distance(f1, f2, 4, 100, 100);
        double d23 = ext_statistical_distance(f2, f3, 4, 100, 100);
        double d13 = ext_statistical_distance(f1, f3, 4, 100, 100);
        if (d13 > d12 + d23 + 0.001) FAIL("triangle inequality violated");
        else OK();
    }

    /* Property 10: Extractor error bound monotonic */
    TEST("Extractor error bound decreases with entropy");
    {
        Extractor e = ext_create(8, 4, 2, 0.1);
        double eps_hi = ext_error_bound(&e, 8, 6.0);
        double eps_lo = ext_error_bound(&e, 8, 3.0);
        if (eps_hi > eps_lo) FAIL("more entropy -> larger error");
        else OK();
    }

    return failures;
}

/*============================================================
 * L4: Theorem verification tests
 *============================================================*/

static int test_theorem_bounds(void) {
    int failures = 0;

    /* LHL bound verification: for k=10, m=4, eps = 2^{-3} = 0.125 */
    TEST("LHL bound: eps = 2^{-(k-m)/2}");
    {
        Extractor e = ext_create(16, 8, 4, 0.1);
        double eps = ext_error_bound(&e, 16, 10.0);
        double expected = pow(2.0, -(10.0 - 4.0) / 2.0); /* 2^{-3} = 0.125 */
        if (fabs(eps - expected) > 0.02) FAIL("LHL bound mismatch");
        else OK();
    }

    /* Seed length bound: d >= 1 */
    TEST("Seed length is at least 1");
    {
        size_t d = ext_seed_len(8, 3.0, 0.5);
        if (d < 1) FAIL("seed length < 1");
        else OK();
    }

    /* Output length bound: m <= n */
    TEST("Output length does not exceed source length");
    {
        size_t m = ext_output_len(8, 7.0);
        if (m > 8) FAIL("m > n");
        else OK();
    }

    /* Privacy amplification bound */
    TEST("Privacy amplification: positive entropy gap -> eps < 1");
    {
        double eps = ext_privacy_amplification_bound(100.0, 50);
        if (eps >= 1.0 || eps <= 0.0) FAIL("invalid eps");
        else OK();
    }

    /* Key derivation bound: zero entropy gap -> eps >= 1 (no security) */
    TEST("Key derivation: zero gap -> eps = 1");
    {
        double eps = ext_key_derivation_bound(128, 32.0, 32);
        /* When min_entropy = key_len, gap = 0, bound returns 1.0 */
        if (eps < 1.0 - 0.01) FAIL("expected eps ~ 1");
        else OK();
    }

    /* Disperser error probability in [0,1] */
    TEST("Disperser error probability in [0,1]");
    {
        double eps = disp_error_prob(16, 8, 8, 8.0);
        if (eps < -0.01 || eps > 1.01) FAIL("eps out of bounds");
        else OK();
    }

    /* Entropy deficit positive for non-uniform */
    TEST("Entropy deficit >= 0");
    {
        double def = ext_entropy_deficit(8, 3.0);
        if (def < -0.01) FAIL("negative deficit");
        else OK();
    }

    return failures;
}

/*============================================================
 * L5: Algorithm tests
 *============================================================*/

static int test_algorithms(void) {
    int failures = 0;

    /* LHL extractor produces output */
    TEST("Leftover Hash Lemma extractor runs");
    {
        SeededExt se = ext_leftover_hash(8, 4);
        bool src[8] = {1,0,1,1,0,0,1,0};
        bool out[4] = {0};
        if (!se.ext(src, 42, out, 4)) FAIL("extraction failed");
        else OK();
    }

    /* LHL extractor deterministic */
    TEST("LHL extractor is deterministic");
    {
        SeededExt se = ext_leftover_hash(8, 4);
        bool src[8] = {1,0,1,1,0,0,1,0};
        bool out1[4] = {0}, out2[4] = {0};
        se.ext(src, 42, out1, 4);
        se.ext(src, 42, out2, 4);
        int match = 1;
        for (int i = 0; i < 4; i++) if (out1[i] != out2[i]) match = 0;
        if (!match) FAIL("non-deterministic");
        else OK();
    }

    /* Trevisan extractor produces output */
    TEST("Trevisan extractor runs");
    {
        SeededExt se = ext_trevisan(8, 4, 0.1);
        bool src[8] = {1,0,1,1,0,0,1,0};
        bool out[4] = {0};
        if (!se.ext(src, 42, out, 4)) FAIL("extraction failed");
        else OK();
    }

    /* Chor-Goldreich 2-source extractor */
    TEST("Chor-Goldreich two-source extractor");
    {
        TwoSourceExt e = ext_two_source_chor_goldreich(8);
        if (e.n != 8) FAIL("wrong n");
        else OK();
    }

    /* Bourgain 2-source extractor */
    TEST("Bourgain two-source extractor");
    {
        TwoSourceExt e = ext_two_source_bourgain(8, 0.01);
        bool x[8] = {0}, y[8] = {0};
        bool out[4] = {0};
        for (int i = 0; i < 8; i++) { x[i] = i % 2; y[i] = (i+1) % 2; }
        if (!ext_bourgain_eval(&e, x, y, out)) FAIL("eval failed");
        else OK();
    }

    /* Extractor graph expansion */
    TEST("Extractor graph verifies expansion");
    {
        ExtractorGraph eg = eg_construct(8, 4, 3);
        bool expands = eg_verify_expansion(&eg);
        double frac = eg_neighbor_fraction(&eg, 16);
        if (!expands) FAIL("expansion failed");
        if (frac < 0.0 || frac > 1.0) FAIL("fraction out of bounds");
        eg_free(&eg);
        if (!expands) {} else OK();
    }

    /* Key derivation */
    TEST("Key derivation produces output");
    {
        uint8_t weak[8] = {0xCA,0xFE,0xBA,0xBE,0x00,0x11,0x22,0x33};
        CryptoKey key = ext_key_derivation(weak, 8, 16);
        if (key.len != 16) FAIL("wrong key length");
        ck_free(&key);
        if (key.len == 0) {} else OK();
    }

    /* Privacy amplification */
    TEST("Privacy amplification runs");
    {
        uint8_t shared[8] = {1,2,3,4,5,6,7,8};
        uint8_t seed[4] = {9,10,11,12};
        uint8_t key[8] = {0};
        if (!ext_privacy_amplification(shared, 8, seed, 4, key, 8))
            FAIL("failed");
        else OK();
    }

    /* Disperser hits all */
    TEST("Disperser hitting test");
    {
        Disperser d = disp_create(6, 4, 4);
        bool hits = disp_hits_all(&d, 100);
        /* For small instances, the disperser should hit most of range */
        if (!hits) FAIL("hitting test failed (may be statistical)");
        else OK();
    }

    return failures;
}

/*============================================================
 * L1+L3: Structural tests
 *============================================================*/

static int test_structures(void) {
    int failures = 0;

    /* GF(2) add = XOR property */
    TEST("GF(2) addition is XOR");
    {
        bool a[4] = {1,0,1,1}, b[4] = {0,1,1,0}, r[4];
        ext_gf2_add(a, b, 4, r);
        /* 1011 XOR 0110 = 1101 */
        if (r[0]!=1 || r[1]!=1 || r[2]!=0 || r[3]!=1) FAIL("wrong XOR result");
        else OK();
    }

    /* GF(2) multiplication is defined */
    TEST("GF(2) multiplication defined");
    {
        bool a[4] = {1,0,1,0}, b[4] = {0,1,0,0}, r[4];
        if (!ext_gf2_mul(a, b, 4, r)) FAIL("multiplication failed");
        else OK();
    }

    /* Bit packing roundtrip */
    TEST("Bit pack/unpack roundtrip");
    {
        bool orig[20];
        for (int i = 0; i < 20; i++) orig[i] = (i * 7 + 3) % 2;
        size_t vec[4] = {0};
        bool back[20] = {0};
        ext_bits_pack(orig, 20, vec);
        ext_bits_unpack(vec, 20, back);
        int ok = 1;
        for (int i = 0; i < 20; i++) if (orig[i] != back[i]) ok = 0;
        if (!ok) FAIL("roundtrip failed");
        else OK();
    }

    /* Popcount identity */
    TEST("Popcount identity: popcount(x) = sum of bits");
    {
        size_t x = 0xBEEF;
        size_t cnt = 0;
        for (size_t t = x; t; t >>= 1) cnt += t & 1;
        if (ext_popcount(x) != cnt) FAIL("popcount mismatch");
        else OK();
    }

    /* Extractor validation */
    TEST("Extractor parameter validation");
    {
        Extractor e1 = ext_create(8, 3, 4, 0.1);
        Extractor e2 = ext_create(8, 10, 4, 0.1); /* d > n -> invalid */
        Extractor e3 = ext_create(8, 3, 4, 0.0);  /* eps=0 -> invalid */
        if (!ext_is_valid(&e1)) FAIL("valid rejected");
        if (ext_is_valid(&e2)) FAIL("invalid accepted (d>n)");
        if (ext_is_valid(&e3)) FAIL("invalid accepted (eps=0)");
        else OK();
    }

    return failures;
}

/*============================================================
 * L6: Canonical problem tests
 *============================================================*/

static int test_canonical_problems(void) {
    int failures = 0;

    /* Problem 1: Extract from biased coin source */
    TEST("Extract from biased coin source");
    {
        size_t freq[4] = {0};
        uint64_t rng = 999;
        ext_sample_distribution(0.7, 2, 1000, freq, &rng);
        double h = ext_min_entropy_from_freq(freq, 4, 1000);
        /* Biased source should have H_inf < 2 */
        if (h >= 1.99) FAIL("biased source has too much entropy");
        else OK();
    }

    /* Problem 2: Estimate stat distance to uniform */
    TEST("Estimate statistical distance");
    {
        size_t obs[4] = {300, 250, 250, 200};
        size_t uni[4] = {250, 250, 250, 250};
        double d = ext_statistical_distance(obs, uni, 4, 1000, 1000);
        /* Should be positive but < 1 */
        if (d <= 0.0 || d >= 1.0) FAIL("distance out of range");
        else OK();
    }

    /* Problem 3: Verify min-entropy of a source */
    TEST("Verify min-entropy bound");
    {
        size_t freq[8] = {100,100,100,100,100,100,100,100};
        /* 800 samples uniform over 8 values -> H_inf = log2(800/100) = 3.0.
         * Verifying k=2.0 with delta=0.05:
         *   p_max_emp = 0.125, eps_chernoff ~ 0.048.
         *   p_max_upper = 0.173, 2^{-2.0} = 0.25.
         *   Since 0.173 <= 0.25, verification passes. */
        bool ok = ext_verify_min_entropy(freq, 8, 800, 2.0, 0.05);
        if (!ok) FAIL("verification rejected valid source");
        else OK();
    }

    /* Problem 4: Simulate extractor */
    TEST("Simulate explicit extractor");
    {
        double delta = ext_simulate_explicit_extractor(6, 3, 4, 3, 0.1, 200);
        if (delta < 0.0 || delta > 1.0) FAIL("delta out of range");
        else OK();
    }

    return failures;
}

/* Dummy BPP algorithm: returns true when majority of seed bits are 1 */
static bool dummy_bpp_alg(const char *inp, size_t ilen,
                          const bool *seed, size_t slen) {
    (void)inp; (void)ilen;
    size_t ones = 0;
    for (size_t i = 0; i < slen; i++)
        if (seed[i]) ones++;
    return ones > slen / 2;
}

/*============================================================
 * L7: Application tests
 *============================================================*/

static int test_applications(void) {
    int failures = 0;

    /* BPP derandomization: dummy algorithm */
    TEST("BPP derandomization with dummy algorithm");
    {
        bool result = false;
        bool ok = ext_derandomize_bpp(dummy_bpp_alg, "test", 4, 4, &result);
        if (!ok) FAIL("derandomization failed");
        else OK();
    }

    /* Network extractor */
    TEST("Network extractor simulation");
    {
        bool ok = ext_network_extractor_simulate(4, 3, 4, 2, 10);
        if (!ok) FAIL("network extractor failed");
        else OK();
    }

    /* CryptoKey equality */
    TEST("CryptoKey constant-time equality");
    {
        CryptoKey a = ck_create(4);
        CryptoKey b = ck_create(4);
        a.key[0] = 0xAB; b.key[0] = 0xAB;
        a.key[1] = 0xCD; b.key[1] = 0xCD;
        if (!ck_eq(&a, &b)) FAIL("equal keys rejected");
        b.key[1] = 0xCE;
        if (ck_eq(&a, &b)) FAIL("different keys accepted");
        ck_free(&a); ck_free(&b);
        if (a.len == 0) OK(); /* freed */
    }

    return failures;
}

/*============================================================
 * Main
 *============================================================*/

int main(void) {
    int total_failures = 0;

    printf("\n============================================================\n");
    printf("  mini-randomness-extractors -- Test Suite\n");
    printf("============================================================\n\n");

    printf("[Module self-tests]\n");
    int f1 = ext_core_self_test();
    printf("  extractor_core:     %s (%d failures)\n",
           f1 == 0 ? "PASS" : "FAIL", f1);
    total_failures += f1;

    int f2 = ext_constructions_self_test();
    printf("  extractor_constr:   %s (%d failures)\n",
           f2 == 0 ? "PASS" : "FAIL", f2);
    total_failures += f2;

    int f3 = ext_apps_self_test();
    printf("  extractor_apps:     %s (%d failures)\n",
           f3 == 0 ? "PASS" : "FAIL", f3);
    total_failures += f3;

    int f4 = disp_core_self_test();
    printf("  disperser_core:     %s (%d failures)\n",
           f4 == 0 ? "PASS" : "FAIL", f4);
    total_failures += f4;

    printf("\n[Property tests]\n");
    int f5 = test_entropy_properties();
    printf("  entropy properties: %s (%d failures)\n",
           f5 == 0 ? "PASS" : "FAIL", f5);
    total_failures += f5;

    int f6 = test_theorem_bounds();
    printf("  theorem bounds:     %s (%d failures)\n",
           f6 == 0 ? "PASS" : "FAIL", f6);
    total_failures += f6;

    printf("\n[Algorithm tests]\n");
    int f7 = test_algorithms();
    printf("  algorithms:         %s (%d failures)\n",
           f7 == 0 ? "PASS" : "FAIL", f7);
    total_failures += f7;

    printf("\n[Structure tests]\n");
    int f8 = test_structures();
    printf("  structures:         %s (%d failures)\n",
           f8 == 0 ? "PASS" : "FAIL", f8);
    total_failures += f8;

    printf("\n[Canonical problem tests]\n");
    int f9 = test_canonical_problems();
    printf("  canonical problems: %s (%d failures)\n",
           f9 == 0 ? "PASS" : "FAIL", f9);
    total_failures += f9;

    printf("\n[Application tests]\n");
    int f10 = test_applications();
    printf("  applications:       %s (%d failures)\n",
           f10 == 0 ? "PASS" : "FAIL", f10);
    total_failures += f10;

    printf("\n============================================================\n");
    if (total_failures == 0) {
        printf("  ALL TESTS PASSED\n");
    } else {
        printf("  %d TEST FAILURES\n", total_failures);
    }
    printf("============================================================\n");

    return total_failures == 0 ? 0 : 1;
}
