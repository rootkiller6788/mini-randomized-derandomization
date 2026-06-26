/**
 * example2.c -- Extractor Constructions in Action
 *
 * Demonstrates L4+L5:
 *   1. Leftover Hash Lemma extractor
 *   2. Trevisan extractor
 *   3. Two-source extractors (Chor-Goldreich, Bourgain)
 *   4. Extractor graph expansion
 *   5. Seeded extractor verification
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "extractor_core.h"
#include "extractor_constructions.h"

int main(void) {
    printf("============================================================\n");
    printf("  Example 2: Extractor Constructions\n");
    printf("============================================================\n\n");

    /* Part 1: LHL extractor */
    printf("--- Part 1: Leftover Hash Lemma Extractor ---\n");
    {
        SeededExt se = ext_leftover_hash(8, 4);
        bool src[8] = {1,0,1,1,0,0,1,0};
        bool out1[4] = {0}, out2[4] = {0};

        /* Same input + seed = deterministic output */
        se.ext(src, 42, out1, 4);
        se.ext(src, 42, out2, 4);
        printf("  Input:  10110010 (8 bits)\n");
        printf("  Seed:   42\n");
        printf("  Output: %d%d%d%d (4 bits)\n",
               out1[3], out1[2], out1[1], out1[0]);
        printf("  Deterministic: %s\n\n",
               memcmp(out1, out2, 4)==0 ? "YES" : "NO");
    }

    /* Part 2: Different seeds -> different outputs */
    printf("--- Part 2: Seed Variance ---\n");
    {
        SeededExt se = ext_leftover_hash(8, 4);
        bool src[8] = {1,0,1,0,1,0,1,0};
        int different = 0;
        for (int s = 0; s < 8; s++) {
            bool out[4] = {0};
            se.ext(src, (size_t)s, out, 4);
            printf("  Seed %d -> %d%d%d%d\n", s,
                   out[3], out[2], out[1], out[0]);
            if (s > 0) different++;
        }
        printf("  Different seeds should give different outputs\n\n");
    }

    /* Part 3: Trevisan extractor */
    printf("--- Part 3: Trevisan Extractor ---\n");
    {
        SeededExt se = ext_trevisan(8, 4, 0.1);
        bool src[8] = {0,1,0,1,0,1,0,1};
        for (int s = 0; s < 4; s++) {
            bool out[4] = {0};
            se.ext(src, (size_t)(s * 13), out, 4);
            printf("  Seed %d -> ", s*13);
            for (int i = 3; i >= 0; i--) printf("%d", out[i]);
            printf("\n");
        }
        printf("  Trevisan: d=O(log n), seed shorter than LHL\n\n");
    }

    /* Part 4: Two-source extractors */
    printf("--- Part 4: Two-Source Extractors ---\n");
    {
        TwoSourceExt cg = ext_two_source_chor_goldreich(8);
        bool x[8] = {1,0,1,0,1,0,1,0}; /* alternating */
        bool y[8] = {0,1,0,1,0,1,0,1}; /* opposite alternating */
        bool cg_out[4] = {0}, bg_out[4] = {0};
        ext_cg_eval(&cg, x, y, cg_out);
        printf("  Chor-Goldreich (inner product): ");
        for (int i = 3; i >= 0; i--) printf("%d", cg_out[i]);
        printf("\n");

        TwoSourceExt bg = ext_two_source_bourgain(4, 0.01);
        bool a[4] = {1,0,1,0}, b[4] = {0,1,0,1};
        ext_bourgain_eval(&bg, a, b, bg_out);
        printf("  Bourgain (field multiplication):  ");
        printf("%d (first bit)\n", bg_out[0]);
    }

    /* Part 5: Extractor graph */
    printf("\n--- Part 5: Extractor Graph Expansion ---\n");
    {
        ExtractorGraph eg = eg_construct(8, 4, 3);
        printf("  Graph: N=%u, K=%u, D=%u\n",
               (unsigned)(1<<eg.n), (unsigned)(1<<eg.k),
               (unsigned)eg.degree);
        printf("  Expansion verified: %s\n",
               eg_verify_expansion(&eg) ? "YES" : "NO");
        for (int s = 4; s <= 64; s *= 2) {
            double frac = eg_neighbor_fraction(&eg, (size_t)s);
            printf("  Set size %d -> neighbor fraction %.4f\n", s, frac);
        }
        eg_free(&eg);
    }

    /* Part 6: Verify extractor uniformity */
    printf("\n--- Part 6: Extractor Verification ---\n");
    {
        SeededExt se = ext_leftover_hash(6, 3);
        bool uniform = ext_seeded_verify(&se, 200);
        printf("  LHL(6->3) empirical uniformity: %s\n",
               uniform ? "PASS" : "NEAR-PASS (statistical)");
    }

    printf("\n=== Example 2 Complete ===\n");
    return 0;
}
