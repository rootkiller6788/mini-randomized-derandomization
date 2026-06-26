/**
 * example2.c — Derandomization using NW PRG
 *
 * Demonstrates:
 *   L7: BPP derandomization via NW PRG
 *   L4: Impagliazzo-Wigderson theorem
 *   L5: Hardness amplification via Yao XOR Lemma
 *
 * This example simulates derandomizing a randomized algorithm
 * for the "String Equality with Noise" problem.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "nw_core.h"
#include "nw_hardness.h"
#include "nw_prg.h"

/* Simulated BPP algorithm: decides if two strings are equal
 * by random sampling of positions.
 * Real BPP algorithm: uses O(log n) random bits for sampling. */
static bool bp_string_eq(const bool *input, const bool *random,
                          size_t n, size_t m) {
    (void)random;
    (void)m;
    /* Trivial: just compare first bits */
    size_t half = n / 2;
    for (size_t i = 0; i < half; i++)
        if (input[i] != input[i + half]) return false;
    return true;
}

int main(void) {
    printf("╔══════════════════════════════════════════╗\n");
    printf("║  BPP Derandomization via NW PRG         ║\n");
    printf("╚══════════════════════════════════════════╝\n\n");

    /* Problem: check if x1 = x2 (two n-bit halves are equal) */
    size_t n = 16;  /* total input size: two 8-bit strings */
    size_t m = 8;   /* "randomness" used by BPP algorithm */

    printf("── Step 1: Hardness Assumption ───────────\n");
    HardnessAssumption ha = {
        .n = 16,
        .circuit_size = 10000.0,
        .hardness = 0.1,
        .is_exponential = false
    };
    printf("Assumption: E requires circuits of size %.0f\n", ha.circuit_size);
    printf("Valid: %s\n\n", nw_hardness_assumption_valid(&ha) ? "YES" : "NO");

    printf("── Step 2: Yao XOR Lemma Analysis ────────\n");
    double eps_base = 0.1;
    for (size_t k = 1; k <= 5; k++) {
        double amplified = nw_yao_xor_amplify(eps_base, k, 0.01);
        printf("  XOR-%zu: epsilon = %.2f → %.4f\n",
               k, eps_base, amplified);
    }
    printf("\n");

    printf("── Step 3: Construct NW PRG for Derandomization ──\n");
    NWPRG *prg = nw_prg_for_BPP(n, m, ha.circuit_size);
    if (!prg) {
        printf("PRG construction failed — insufficient hardness.\n");
        printf("This demonstrates the hardness assumption requirement.\n");
        return 0;
    }

    printf("PRG constructed: seed_len=%zu, output_len=%zu\n",
           prg->seed_len, prg->output_len);
    nw_prg_print(prg, stdout);
    printf("\n");

    printf("── Step 4: Derandomize BPP Algorithm ─────\n");
    /* Input: two equal strings */
    bool input[16];
    for (size_t i = 0; i < 8; i++) {
        input[i] = (i % 2) ? true : false;
        input[i + 8] = (i % 2) ? true : false;  /* equal */
    }

    bool result;
    bool ok = nw_bpp_derandomize(prg, bp_string_eq, input, n, m, &result);
    printf("BPP algorithm deterministically simulated.\n");
    printf("Result: %s (expected: EQUAL)\n",
           result ? "EQUAL" : "NOT EQUAL");
    printf("Correct: %s\n\n", result ? "YES" : "NO");

    printf("── Step 5: IW97 Derandomization Levels ───\n");
    const char *levels[] = {"SUBEXP", "QP", "P"};
    for (int lvl = 0; lvl < 3; lvl++) {
        size_t seed;
        bool feasible = nw_iw_derandomize(lvl, n, &seed);
        printf("Level %d (%s): feasible=%s, seed_len=%zu\n",
               lvl, levels[lvl],
               feasible ? "YES" : "NO", seed);
    }

    printf("\n── Done ──────────────────────────────────\n");

    if (prg) nw_prg_free(prg);
    return 0;
}
