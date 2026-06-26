/*
 * example3.c — Goldreich-Levin Algorithm and Hardness Amplification
 *
 * Demonstrates:
 *   1. Creating a Hadamard code (truth table of linear functions)
 *   2. Encoding a message as a Hadamard codeword
 *   3. Goldreich-Levin list-decoding with oracle
 *   4. Fourier analysis of Boolean functions
 *   5. Hard-core predicate extraction
 *   6. Hardness amplification via list-decodable codes
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "hadamard.h"
#include "reed_muller.h"
#include "list_decode_apps.h"

/* Simulated oracle: agrees with linear function <x, y> = <101, y> on ~80% of inputs */
static bool noisy_oracle(const bool *y) {
    bool inner = (y[0] && true) ^ (y[1] && false) ^ (y[2] && true);
    /* Add noise: flip with probability ~0.15 */
    /* Use a simple pseudo-random test based on input bits */
    bool noise = (y[0] != y[1]) && !y[2];
    return inner ^ noise;
}

static bool pure_oracle(const bool *y) {
    return (y[0] && true) ^ (y[1] && false) ^ (y[2] && true);
}

int main(void) {
    printf("=== Goldreich-Levin Algorithm & Hardness Amplification ===\n\n");

    /* Step 1: Hadamard code H_3 */
    printf("Step 1: Creating Hadamard code H_3 (k=3, n=8, d=4)...\n");
    HadamardCode *hc = had_create(3);
    if (!hc) {
        printf("  Failed.\n");
        return 1;
    }
    printf("  n=%zu, k=%zu, d=%zu (rate=%.4f)\n",
           hc->n, hc->k, hc->d, (double)hc->k / (double)hc->n);

    /* Step 2: Encode message x = (1, 0, 1) */
    printf("\nStep 2: Encoding message x = (1, 0, 1)...\n");
    bool msg[] = {true, false, true};
    bool codeword[8];
    had_encode(hc, msg, codeword);
    printf("  Hadamard codeword: [");
    for (size_t i = 0; i < 8; i++) printf("%d ", codeword[i] ? 1 : 0);
    printf("]\n");

    /* Step 3: Goldreich-Levin list-decoding */
    printf("\nStep 3: Goldreich-Levin algorithm (epsilon=0.25)...\n");
    bool *gl_out = NULL;
    size_t L = had_goldreich_levin(hc, noisy_oracle, 0.25, &gl_out, 10);
    printf("  Found %zu candidate(s)\n", L);
    if (L > 0 && gl_out) {
        for (size_t l = 0; l < L; l++) {
            printf("  Candidate %zu: [%d, %d, %d]\n", l,
                   gl_out[l * hc->k + 0] ? 1 : 0,
                   gl_out[l * hc->k + 1] ? 1 : 0,
                   gl_out[l * hc->k + 2] ? 1 : 0);
        }
        free(gl_out);
    }

    /* Step 4: Fourier analysis */
    printf("\nStep 4: Fourier analysis of the oracle...\n");
    double coeff = had_fourier_coefficient(pure_oracle, 3, msg);
    printf("  Fourier coefficient at (1,0,1): %.4f\n", coeff);
    double abs_coeff = fabs(coeff);
    printf("  |f_hat| >= 0.5: %s (heavy Fourier coefficient)\n",
           abs_coeff >= 0.5 ? "YES" : "NO");

    /* Step 5: Hard-core predicate */
    printf("\nStep 5: Hard-core predicate extraction (Goldreich-Levin)...\n");
    bool r[] = {true, true, false};
    bool hcb = had_hard_core_bit(msg, r, 3);
    printf("  GL(x=%d%d%d, r=%d%d%d) = %d\n",
           msg[0]?1:0, msg[1]?1:0, msg[2]?1:0,
           r[0]?1:0, r[1]?1:0, r[2]?1:0,
           hcb ? 1 : 0);

    /* Step 6: Hardness amplification via codes */
    printf("\nStep 6: Hardness amplification via list-decodable codes...\n");
    HardnessAmpParams hap;
    hap.n = 8;
    hap.eps = 0.15;
    hap.L = 4;
    double amplified = 0.0;
    bool dummy_fn[256] = {0};
    ld_hardness_amplification(dummy_fn, hap.n, &hap, &amplified);
    printf("  Original (1-eps)-hard, eps=%.2f\n", hap.eps);
    printf("  Amplified hardness: %.4f (target: > 0.5+threshold)\n",
           amplified);

    /* Johnson bound for list decoding */
    printf("\nStep 7: List-size bounds...\n");
    double eps = 0.10;
    size_t ls_bound = had_list_size_bound(eps);
    printf("  Hadamard list-size bound (eps=%.2f): L <= %zu\n",
           eps, ls_bound);

    printf("\n=== Example Complete ===\n");

    had_free(hc);
    return 0;
}
