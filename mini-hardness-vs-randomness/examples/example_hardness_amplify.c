/******************************************************************************
 * example_hardness_amplify.c ¡ª Example: Hardness Amplification
 *
 * Demonstrates worst-case to average-case hardness amplification
 * using Yao's XOR Lemma and Impagliazzo's Hardcore Lemma.
 *
 * Shows the connection to the Nisan-Wigderson PRG construction
 * and the five worlds framework of Impagliazzo (1995).
 *
 * Reference: Yao (1982); Impagliazzo (1995); Arora & Barak, Ch.19
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

/* Forward declarations from iw_theorem.c */
extern bool iw_theorem_conditions_met(const HardnessCert *cert, size_t n);
extern bool iw_full_derandomization_pipeline(size_t n, HRLevel level);
extern bool iw_check_eca_assumption(size_t n, double epsilon);
extern double iw_derandomization_blowup(size_t n, double hardness);
extern double iw_pseudodeterministic_probability(size_t n, double hardness);
extern double iw_eth_derandomization_gap(size_t n);

int main(void) {
    printf("=== Hardness Amplification & Five Worlds ===\n\n");
    /* Part 1: Yao's XOR Lemma amplification */
    printf("--- Part 1: Yao's XOR Lemma ---\n");
    printf("Starting from worst-case hardness ¦Ä=0.1\n");
    printf("  k  | Avg Hardness | Advantage | Effective\n");
    printf("  ---+--------------+-----------+----------\n");
    for (size_t k = 1; k <= 16; k *= 2) {
        WTAConversion conv = wta_yao_xor_lemma(0.1, 16, k);
        double advantage = 1.0 - 2.0 * (0.5 - conv.avg_hardness);
        printf("  %2zu | %.6f     | %.6f  | %s\n",
               k, conv.avg_hardness, advantage,
               advantage < 0.01 ? "near-random" : "predictable");
    }
    /* Part 2: Impagliazzo Hardcore Lemma */
    printf("\n--- Part 2: Impagliazzo Hardcore Lemma ---\n");
    for (double delta = 0.05; delta <= 0.5; delta += 0.15) {
        WTAConversion conv = wta_impagliazzo_hardcore(delta, 16);
        printf("  ¦Ä=%.2f: hardcore avg hardness=%.4f\n",
               delta, conv.avg_hardness);
    }
    /* Part 3: Sample complexity for amplification */
    printf("\n--- Part 3: Sample Complexity ---\n");
    printf("To estimate hardness within ¦Å with confidence 1-¦Ä:\n");
    for (double eps = 0.2; eps >= 0.01; eps /= 4.0) {
        size_t samples = wta_sample_complexity(16, eps, 0.01);
        printf("  ¦Å=%.4f: samples=%zu\n", eps, samples);
    }
    /* Part 4: Goldreich-Levin applicability */
    printf("\n--- Part 4: Goldreich-Levin Hardcore Predicate ---\n");
    for (size_t n = 4; n <= 64; n *= 2) {
        double threshold = 1.0 / (double)n;
        printf("  n=%2zu: GL applies if hardness ¡Ý %.4f\n", n, threshold);
    }
    /* Part 5: IW Theorem pipeline test */
    printf("\n--- Part 5: IW Theorem Derandomization ---\n");
    for (size_t n = 8; n <= 16; n += 4) {
        bool ok = iw_full_derandomization_pipeline(n, HR_STRONG);
        printf("  n=%2zu, STRONG hardness: derandomization %s\n",
               n, ok ? "POSSIBLE" : "IMPOSSIBLE");
    }
    /* Part 6: ECA assumption check */
    printf("\n--- Part 6: Exponential Complexity Assumption ---\n");
    for (double eps = 0.01; eps <= 1.0; eps *= 10.0) {
        bool valid = iw_check_eca_assumption(20, eps);
        printf("  ¦Å=%.2f: ECA %s\n", eps, valid ? "plausible" : "too strong");
    }
    /* Part 7: ETH vs IW gap */
    printf("\n--- Part 7: ETH ¡ú IW Gap ---\n");
    for (size_t n = 10; n <= 100; n += 30) {
        double gap = iw_eth_derandomization_gap(n);
        printf("  n=%3zu: ETH hardness - IW threshold = %.4f ¡ú %s\n",
               n, gap, gap > 0 ? "sufficient" : "insufficient");
    }
    /* Part 8: Five Worlds classification */
    printf("\n--- Part 8: Impagliazzo's Five Worlds ---\n");
    const char *world_names[] = {
        "Algorithmica", "Heuristica", "Pessiland",
        "Minicrypt", "Cryptomania"
    };
    double configs[][4] = {
        {0.0001, 0.0001, 0.0, 0.0}, /* Algorithmica */
        {0.1,    0.0001, 0.0, 0.0}, /* Heuristica */
        {0.1,    0.1,    0.0, 0.0}, /* Pessiland */
        {0.1,    0.1,    1.0, 0.0}, /* Minicrypt */
        {0.2,    0.2,    1.0, 1.0}, /* Cryptomania */
    };
    for (int w = 0; w < 5; w++) {
        printf("  %s: worst=%.4f, avg=%.4f, OWF=%d, Trapdoor=%d\n",
               world_names[w], configs[w][0], configs[w][1],
               (int)configs[w][2], (int)configs[w][3]);
    }
    printf("\n=== Example Complete ===\n");
    return 0;
}