/**
 * example2.c - Nisan's PRG: Derandomizing Space-Bounded Computation
 * ==================================================================
 * Demonstrates Nisan's pseudorandom generator for space-bounded
 * computation and BPL derandomization.
 *
 * Knowledge: L5 (Nisan PRG algorithm), L7 (BPL derandomization).
 * Refs: Nisan (1992), Arora-Barak Section 20.3.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "space_derand.h"
#include "nisan_prg.h"

/* Simple test function: counts parity of first 8 bits */
static bool parity_test(const bool *bits, size_t len, void *ctx) {
    (void)ctx;
    size_t limit = len < 8 ? len : 8;
    bool parity = false;
    for (size_t i = 0; i < limit; i++) parity ^= bits[i];
    return parity;
}

int main(void) {
    printf("=== Nisan's PRG: Derandomization Example ===\n\n");

    /* PRG Construction */
    printf("--- PRG Construction ---\n");
    size_t space = 8;       /* O(log n) space */
    size_t runtime = 256;   /* poly(n) time steps */
    size_t seed_len = nisan_prg_seed_len(space, runtime);
    size_t output_len = nisan_prg_output_len(seed_len, space);
    printf("  Space S=%zu, Time T=%zu\n", space, runtime);
    printf("  Seed length: %zu bits\n", seed_len);
    printf("  Output length: %zu bits\n", output_len);
    printf("  Stretch factor: %.2fx\n",
           (double)output_len / (double)seed_len);

    /* Create and evaluate PRG */
    printf("\n--- PRG Evaluation ---\n");
    SpacePRG *prg = nisan_prg_create(space, 128);
    if (prg) {
        bool seed[64] = {0};
        bool out[256] = {0};
        /* Set a non-zero seed */
        seed[0] = true; seed[3] = true; seed[7] = true;
        if (nisan_prg_eval(prg, seed, prg->seed_len, out, prg->output_len)) {
            printf("  First 16 output bits: ");
            for (size_t i = 0; i < 16 && i < prg->output_len; i++)
                printf("%d", out[i] ? 1 : 0);
            printf("\n");
        }

        /* Uniformity check */
        printf("\n--- Uniformity Check ---\n");
        bool uniform = nisan_prg_verify_uniformity(prg, 1000);
        printf("  PRG passes uniformity: %s\n", uniform ? "yes" : "no");
        nisan_prg_free(prg);
    }

    /* Hash Family */
    printf("\n--- Pairwise Independent Hash Family ---\n");
    PairwiseHashFamily *phf = nisan_hash_create(8, 8);
    if (phf) {
        printf("  Family size: %zu functions\n", phf->num_functions);
        size_t h0 = nisan_hash_eval(phf, 0, 42);
        size_t h1 = nisan_hash_eval(phf, 1, 42);
        printf("  h_0(42) = %zu, h_1(42) = %zu\n", h0, h1);
        bool pi = nisan_hash_pairwise_independent(phf, 500);
        printf("  Pairwise independent: %s\n", pi ? "yes" : "no");
        nisan_hash_free(phf);
    }

    /* BPL Derandomization */
    printf("\n--- BPL Derandomization ---\n");
    bool result = false;
    bool ok = nisan_derandomize_bpl(space, runtime, parity_test, NULL, &result);
    printf("  Derandomization succeeded: %s\n", ok ? "yes" : "no");
    printf("  Result: %s\n", result ? "accept" : "reject");

    /* Trade-off analysis */
    printf("\n--- Seed-Error Trade-off ---\n");
    for (size_t sl = 10; sl <= 50; sl += 10) {
        double err = nisan_error_for_seed_len(space, 128, sl);
        printf("  seed=%zu: error <= %.6f\n", sl, err);
    }

    printf("\n=== Example2 Complete ===\n");
    return 0;
}
