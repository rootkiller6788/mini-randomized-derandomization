/******************************************************************************
 * example_circuit_lb.c ¡ª Example: Circuit Lower Bound Computation
 *
 * Demonstrates the computation and comparison of various circuit lower
 * bounds: Shannon counting, H?stad switching lemma, Razborov monotone,
 * and Smolensky algebraic bounds.
 *
 * Shows the hardness levels derived from each bound and their
 * implications for derandomization.
 *
 * Reference: Arora & Barak (2009), Chapters 6, 14
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "circuit_lower_bounds.h"
#include "hardness_randomness.h"

int main(void) {
    printf("=== Circuit Lower Bounds for Hardness-vs-Randomness ===\n\n");
    printf("This example computes circuit lower bounds for various\n");
    printf("Boolean functions and circuit classes, demonstrating\n");
    printf("the hardness required for derandomization.\n\n");
    printf("--- L4: Shannon's Counting Lower Bound (1949) ---\n");
    for (size_t n = 2; n <= 12; n += 2) {
        CircuitLB lb = clb_shannon_counting(n);
        double quality = clb_quality_metric(&lb);
        printf("  n=%2zu: lb=%.2e gates, quality=%.3f, type=%s\n",
               n, lb.lower_bound, quality, clb_type_name(&lb));
    }
    printf("\n--- L4: H?stad's Switching Lemma (1986) ---\n");
    printf("  PARITY requires exponential-size AC0 circuits\n");
    for (size_t depth = 2; depth <= 5; depth++) {
        CircuitLB lb = clb_hastad_switching(16, depth);
        printf("  depth=%zu: lb¡Ö%.1f gates\n", depth, lb.lower_bound);
    }
    printf("\n--- L4: Razborov's Monotone Circuit Bound (1985) ---\n");
    printf("  CLIQUE requires super-polynomial monotone circuits\n");
    for (size_t n = 10; n <= 50; n += 10) {
        CircuitLB lb = clb_razborov_clique(n);
        printf("  n=%2zu: lb¡Ö%.1f gates\n", n, lb.lower_bound);
    }
    printf("\n--- L4: Smolensky's Algebraic Bound (1987) ---\n");
    printf("  MODp requires exponential AC0[q] circuits for p¡Ùq\n");
    CircuitLB lb3 = clb_smolensky_modp(32, 3);
    CircuitLB lb5 = clb_smolensky_modp(32, 5);
    printf("  p=3: %.1f gates, p=5: %.1f gates\n", lb3.lower_bound, lb5.lower_bound);
    printf("\n--- Comparing Hardness Levels ---\n");
    CircuitLB sh = clb_shannon_counting(10);
    printf("  Shannon (n=10): STRONG enough? %s\n",
           clb_is_strong_enough(&sh, HR_STRONG) ? "YES" : "NO");
    printf("  Quality metric: %.4f\n", clb_quality_metric(&sh));
    double hardness = hr_hardness_from_circuit_lb(10, (size_t)sh.lower_bound);
    printf("  Hardness parameter: %.4f\n", hardness);
    printf("  Amplification possible: %s\n",
           hr_hardness_amplification_possible(hardness, 10) ? "YES" : "NO");
    printf("\n--- Hardness ¡ú BPP=P Check ---\n");
    const char *levels[] = {"WEAK","MODERATE","STRONG","EXPONENTIAL"};
    HRLevel lvls[] = {HR_WEAK, HR_MODERATE, HR_STRONG, HR_EXPONENTIAL};
    for (int i = 0; i < 4; i++) {
        printf("  %s ¡ú BPP=P: %s\n", levels[i],
               hr_bpp_equals_p_check(lvls[i], 50) ? "YES" : "NO");
    }
    printf("\n=== Example Complete ===\n");
    return 0;
}