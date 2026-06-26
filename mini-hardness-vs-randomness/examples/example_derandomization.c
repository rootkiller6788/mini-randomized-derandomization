/******************************************************************************
 * example_derandomization.c ¡ª Example: Derandomization Pipeline
 *
 * Demonstrates the full derandomization pipeline:
 *   1. Circuit lower bound ¡ú Hardness certificate
 *   2. Hardness certificate ¡ú PRG parameters
 *   3. PRG-based derandomization ¡ú Deterministic result
 *
 * Includes a simulated BPP algorithm (primality testing style)
 * and shows the derandomization in action.
 *
 * Reference: Impagliazzo & Wigderson (1997); Arora & Barak, Ch.20
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "circuit_lower_bounds.h"
#include "hardness_randomness.h"
#include "derandomization_via_hardness.h"
#include "worst_to_average.h"

/* Simulated BPP algorithm: Majority-of-three voting with random bits.
   This mimics a BPP algorithm that uses 3 random bits to decide.
   The "actual" answer is hardcoded for demonstration. */
bool simulated_bpp_alg(const char *input, size_t ilen,
                        const bool *random_bits, size_t rlen) {
    (void)input; (void)ilen;
    /* Decision: if first 3 random bits have majority 1, return 1 */
    if (rlen < 3) return false;
    size_t sum = (size_t)random_bits[0] + (size_t)random_bits[1] + (size_t)random_bits[2];
    return sum >= 2;
}

int main(void) {
    printf("=== Derandomization via Hardness: Full Pipeline ===\n\n");
    printf("This example demonstrates:\n");
    printf("  1) Computing circuit lower bounds\n");
    printf("  2) Creating hardness certificates\n");
    printf("  3) Performing derandomization\n\n");
    /* Step 1: Obtain a hardness certificate */
    size_t n = 10;
    CircuitLB lb = clb_shannon_counting(n);
    printf("Step 1: Circuit lower bound for n=%zu\n", n);
    printf("  Theorem: %s\n", lb.theorem);
    printf("  Lower bound: %.1f gates\n", lb.lower_bound);
    double quality = clb_quality_metric(&lb);
    printf("  Quality: %.4f (log2(lb)/n)\n", quality);
    /* Step 2: Create hardness certificate */
    HardnessCert cert = derand_cert_from_circuit_lb(n, (size_t)lb.lower_bound);
    printf("\nStep 2: Hardness certificate\n");
    printf("  Hardness = log2(%.0f)/%zu = %.4f\n",
           lb.lower_bound, n, cert.hardness);
    printf("  Verified: %s\n",
           derand_verify_certificate(&cert, n) ? "YES" : "NO");
    /* Step 3: Check PRG parameters */
    printf("\nStep 3: PRG parameters (NW construction)\n");
    for (int level = HR_WEAK; level <= HR_EXPONENTIAL; level++) {
        size_t seed = hr_seed_len(50, (HRLevel)level);
        size_t output = hr_output_len(50, (HRLevel)level);
        double stretch = (seed > 0) ? (double)output / (double)seed : 0.0;
        printf("  %s: seed=%zu, output=%zu, stretch=%.2f\n",
               hr_level_name((HRLevel)level), seed, output, stretch);
    }
    /* Step 4: Time complexity of derandomization */
    printf("\nStep 4: Derandomization time blowup\n");
    for (double h = 0.05; h <= 1.0; h *= 2.0) {
        size_t blowup = derand_time_complexity(100, h);
        printf("  hardness=%.2f: blowup=2^%.0f ¡Ö %.1e\n",
               h, log2((double)blowup), (double)blowup);
    }
    /* Step 5: Error probability analysis */
    printf("\nStep 5: Error probability vs trials\n");
    for (size_t t = 10; t <= 100; t += 30) {
        double err = derand_error_prob(n, t, 0.3);
        printf("  trials=%3zu: error=%.6f\n", t, err);
    }
    /* Step 6: Advice size for P/poly */
    printf("\nStep 6: Adleman advice size (BPP ? P/poly)\n");
    for (size_t nn = 8; nn <= 64; nn *= 2) {
        size_t adv = derand_advice_size(nn, 0.5);
        printf("  n=%2zu: advice=%zu bits\n", nn, adv);
    }
    /* Step 7: Worst-to-average conversion */
    printf("\nStep 7: Hardness amplification (Yao XOR Lemma)\n");
    for (size_t k = 1; k <= 8; k *= 2) {
        WTAConversion conv = wta_yao_xor_lemma(0.1, 16, k);
        printf("  k=%zu copies: avg_hardness=%.4f\n", k, conv.avg_hardness);
    }
    printf("\n=== Example Complete ===\n");
    return 0;
}