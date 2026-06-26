/**
 * example_nw_demo.c - Complete NW PRG Pipeline Demo
 *
 * Demonstrates the full Nisan-Wigderson pipeline.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "nw_core.h"
#include "nw_designs.h"
#include "nw_hardness.h"
#include "nw_prg.h"
#include "nw_circuits.h"

int main(void) {
    printf("=== Complete NW PRG Pipeline Demo ===\n\n");
    printf("Nisan & Wigderson (1994): Hardness => Randomness\n");
    printf("Impagliazzo & Wigderson (1997): BPP derandomization\n\n");

    printf("-- Phase 1: Hard Function --\n");
    TruthTable *tt = nw_truth_table_create(nw_default_hard_function, 4);
    assert(tt != NULL);
    printf("Truth table for IP-4: %zu entries\n", tt->size);
    printf("Hardness check: %s\n",
           nw_truth_table_is_hard(tt, 100.0) ? "PASS" : "FAIL");
    nw_truth_table_free(tt);

    printf("\n-- Phase 2: Design Constructions --\n");
    NWDesign *rs = nw_reed_solomon_design(4, 1);
    if (rs) {
        printf("Reed-Solomon: k=%zu, m=%zu, l=%zu\n", rs->k, rs->m, rs->l);
        printf("  Verified: %s\n", nw_design_verify(rs) ? "YES" : "NO");
        nw_design_free(rs);
    }
    NWDesign *rd = nw_random_design(25, 8, 5);
    if (rd) {
        printf("Random design: k=%zu, m=%zu, l=%zu\n", rd->k, rd->m, rd->l);
        nw_design_free(rd);
    }

    printf("\n-- Phase 3: PRG Construction --\n");
    const char *mode_names[] = {"STANDARD", "YAO_XOR", "IW97", "OPTIMIZED"};
    for (int mode = 0; mode < 4; mode++) {
        NWPRG *prg = nw_prg_construct_optimal(64, 1000.0, (NWPRGMode)mode);
        if (prg) {
            printf("%s: seed=%zu, output=%zu, stretch=%.2f\n",
                   mode_names[mode], prg->seed_len, prg->output_len, prg->stretch);
            nw_prg_free(prg);
        }
    }

    printf("\n-- Phase 4: Statistical Testing --\n");
    NWPRG *prg_test = nw_prg_construct_optimal(128, 10000.0, NW_PRG_STANDARD);
    if (prg_test) {
        StatisticalTestResult *stats = nw_prg_statistical_tests(prg_test, 10, 256);
        if (stats) {
            printf("monobit=%.4f, runs=%.4f, longest_run=%.4f\n",
                   stats->p_values[0], stats->p_values[1], stats->p_values[2]);
            printf("Passed: %s\n", stats->passes_all ? "YES" : "NO");
            nw_statistical_test_result_free(stats);
        }
        nw_prg_free(prg_test);
    }

    printf("\n-- Phase 5: BPP Simulation --\n");
    for (double eps = 0.1; eps >= 0.001; eps /= 10.0) {
        size_t samples = nw_bpp_simulation_samples(100, eps);
        printf("  Error %.04f -> %zu samples\n", eps, samples);
    }

    printf("\n-- Phase 6: Circuit Lower Bounds --\n");
    for (size_t n = 2; n <= 8; n *= 2) {
        double lb = nw_shannon_lower_bound(n);
        double ub = nw_lupanov_upper_bound(n);
        printf("  n=%zu: Shannon LB=%.1f, Lupanov UB=%.1f\n", n, lb, ub);
    }

    printf("\n=== Pipeline Complete ===\n");
    return 0;
}
