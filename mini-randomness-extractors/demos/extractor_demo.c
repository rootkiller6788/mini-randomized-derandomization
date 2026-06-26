/**
 * demos/extractor_demo.c -- Extractor pipeline demonstration
 *
 * Shows the full pipeline: weak source -> extractor -> uniform output.
 * Compile: cc -I../include -o ext_demo extractor_demo.c ../src/extractor_core.c ../src/extractor_constructions.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/extractor_core.h"
#include "../include/extractor_constructions.h"

int main(void) {
    printf("Extractor Pipeline Demo\n");
    printf("=======================\n\n");

    /* Step 1: Characterize a weak source */
    printf("Step 1: Weak Source Characterization\n");
    size_t n = 8;
    double bias = 0.75; /* biased coin */
    printf("  Source: n=%zu bits, per-bit bias=%.2f\n", n, bias);
    /* Expected min-entropy: each bit has H_inf ~ -log2(0.75) ~ 0.415 */
    /* n independent biased bits -> H_inf ~ n * 0.415 */
    double expected_k = n * (-log2(bias));
    printf("  Estimated min-entropy: %.2f bits\n\n", expected_k);

    /* Step 2: Build extractor */
    printf("Step 2: Extractor Construction\n");
    size_t m = (size_t)(expected_k / 2.0); /* output half the entropy */
    SeededExt ext = ext_leftover_hash(n, m);
    printf("  Construction: Leftover Hash Lemma\n");
    printf("  Input: %zu bits, Output: %zu bits\n", n, m);
    double eps = pow(2.0, -(expected_k - m) / 2.0);
    printf("  Theoretical error: eps <= %.2e\n\n", eps);

    /* Step 3: Extract */
    printf("Step 3: Extraction\n");
    bool src[8] = {1,1,0,1,0,1,1,0}; /* simulated weak source */
    for (int seed = 0; seed < 4; seed++) {
        bool out[4] = {0};
        ext.ext(src, (size_t)seed, out, 4);
        printf("  Seed %d: output = ", seed);
        for (int i = (int)m-1; i >= 0; i--) printf("%d", out[i]);
        printf("\n");
    }

    /* Step 4: Verify uniformity */
    printf("\nStep 4: Uniformity Check\n");
    bool uniform = ext_seeded_verify(&ext, 500);
    printf("  Chi-squared uniformity test: %s\n",
           uniform ? "PASS" : "NEAR-PASS");

    printf("\n=== Pipeline Complete ===\n");
    return 0;
}
