/******************************************************************************
 * demo.c ¡ª Interactive Hardness vs Randomness Demonstration
 *
 * Interactive demo exploring the hardness-randomness tradeoff.
 * Shows how different hardness levels affect derandomization
 * feasibility and PRG parameters.
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

int main(void) {
    printf("\n");
    printf("  ¨X¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨[\n");
    printf("  ¨U   Hardness vs Randomness ¡ª Interactive Demo     ¨U\n");
    printf("  ¨U   Nisan-Wigderson / Impagliazzo-Wigderson      ¨U\n");
    printf("  ¨^¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨T¨a\n\n");
    printf("Core Theorem: If hard functions exist, then\n");
    printf("pseudorandom generators exist, which can derandomize\n");
    printf("all of BPP (P = BPP).\n\n");
    printf("--- Circuit Lower Bounds ---\n");
    for (size_t n = 4; n <= 20; n += 4) {
        CircuitLB lb = clb_shannon_counting(n);
        printf("  n=%2zu | Shannon lb = %.1e | H?stad = %.1e | Razborov = %.1e\n",
               n, lb.lower_bound,
               clb_hastad_switching(n, 4).lower_bound,
               clb_razborov_clique(n).lower_bound);
    }
    printf("\n--- Hardness Levels ¡ú PRG Stretch ---\n");
    printf("  Level         | Seed    | Output  | BPP=P?\n");
    printf("  --------------+---------+---------+-------\n");
    for (int level = 0; level < 4; level++) {
        HRLevel lv = (HRLevel)level;
        size_t seed = hr_seed_len(50, lv);
        size_t out = hr_output_len(50, lv);
        printf("  %-14s | %7zu | %7zu | %s\n",
               hr_level_name(lv), seed, out,
               hr_bpp_equals_p_check(lv, 50) ? "YES" : "NO");
    }
    printf("\n--- Yao's XOR Lemma: Hardness Amplification ---\n");
    printf("  Starting from ¦Ä=0.1 worst-case hardness:\n");
    for (size_t k = 1; k <= 16; k *= 2) {
        double bound = wta_xor_lemma_bound(0.1, k);
        printf("  k=%2zu ¡ú prediction prob ¡Ü %.4f\n", k, bound);
    }
    printf("\n--- IW Theorem: BPP Derandomization ---\n");
    printf("  Under exponential hardness (E requires 2^{¦¸(n)} circuits):\n");
    printf("  ¡ú BPP = P (all randomized poly-time = deterministic poly-time)\n");
    printf("  ¡ú PRG exists with seed O(log n), output poly(n)\n");
    printf("  ¡ú Cryptographic PRGs from any one-way function\n\n");
    return 0;
}