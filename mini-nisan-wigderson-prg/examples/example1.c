/**
 * example1.c — NW PRG Construction and Evaluation Demo
 *
 * Demonstrates:
 *   L1: Core type definitions
 *   L2: PRG construction from hard function + design
 *   L5: NW PRG evaluation
 *
 * This example constructs a small NW PRG and evaluates it
 * on multiple seeds, showing the pseudorandom output.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nw_core.h"
#include "nw_designs.h"
#include "nw_hardness.h"
#include "nw_prg.h"

/* Simple hard function: PARITY of input bits */
static bool parity_fn(const bool *input, size_t n) {
    bool result = false;
    for (size_t i = 0; i < n; i++)
        result ^= input[i];
    return result;
}

int main(void) {
    printf("╔══════════════════════════════════════════╗\n");
    printf("║  NW PRG Construction & Evaluation Demo  ║\n");
    printf("╚══════════════════════════════════════════╝\n\n");

    /* Step 1: Create a combinatorial design.
     * A (16, 8, 4)-design: 8 subsets of [16], each of size 4,
     * with pairwise intersection at most 3. */
    printf("── Step 1: Combinatorial Design Creation ──\n");
    size_t k = 16, m = 8, l = 4, bound = 3;
    NWDesign *design = nw_design_create(k, m, l, bound);
    assert(design != NULL);

    /* Manually set design sets (disjoint for simplicity):
     * S_0={0,1,2,3}, S_1={4,5,6,7}, ..., S_7={12,13,14,15} */
    for (size_t i = 0; i < m; i++)
        for (size_t j = 0; j < l; j++)
            design->sets[i][j] = i * l + j;

    printf("Design parameters: k=%zu, m=%zu, l=%zu\n", k, m, l);
    printf("Intersection bound: %zu\n", bound);
    printf("Verification: %s\n\n",
           nw_design_verify(design) ? "PASS" : "FAIL");

    /* Step 2: Create hard function.
     * PARITY requires exponential-size circuits in AC0,
     * making it a standard candidate hard function. */
    printf("── Step 2: Hard Function ──────────────────\n");
    HardFunction *hf = nw_hard_function_create(
        4, parity_fn, 100.0, 0.1, "PARITY-4"
    );
    assert(hf != NULL);
    printf("Hard function: %s\n", hf->description);
    printf("Input length: n=%zu, hardness=%.0f\n\n", hf->n, hf->hardness);

    /* Step 3: Construct NW PRG.
     * G: {0,1}^16 → {0,1}^8
     * G(x)_i = f(x|_{S_i}) */
    printf("── Step 3: NW PRG Construction ────────────\n");
    NWPRG *prg = nw_prg_construct(design, hf);
    assert(prg != NULL);
    printf("PRG: {0,1}^%zu → {0,1}^%zu (stretch=%.2f)\n\n",
           prg->seed_len, prg->output_len, prg->stretch);

    /* Step 4: Evaluate PRG on sample seeds */
    printf("── Step 4: PRG Evaluation ────────────────\n");
    printf("Seed (hex)  →  Output (hex)  |  Notes\n");
    printf("──────────     ─────────────  ────────────\n");

    uint64_t test_seeds[] = {0x0000, 0x0001, 0x00FF, 0x5555, 0xAAAA};
    for (int s = 0; s < 5; s++) {
        uint64_t seed_val = test_seeds[s];
        bool seed[16];
        nw_int_to_seed(16, seed_val, seed);

        bool output[8];
        nw_prg_evaluate(prg, seed, output);

        uint64_t out_val = 0;
        for (size_t i = 0; i < 8; i++)
            if (output[i]) out_val |= (1ULL << i);

        printf("  0x%04zX      →  0x%02zX          |  ", seed_val, out_val);
        if (seed_val == 0)
            printf("zero seed\n");
        else if (seed_val == 0x5555)
            printf("alternating pattern\n");
        else
            printf("\n");
    }

    printf("\n── Step 5: Seed Length Analysis ───────────\n");
    for (size_t out_len = 10; out_len <= 1000; out_len *= 10) {
        size_t seed_len = nw_seed_length(out_len, 1000.0);
        printf("Output %4zu bits → Seed %3zu bits (ratio %.1fx)\n",
               out_len, seed_len, (double)out_len / (double)seed_len);
    }

    printf("\n── Done ──────────────────────────────────\n");
    nw_prg_free(prg);
    return 0;
}
