#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "prg_core.h"
#include "prg_constructions.h"
#include "hardcore_predicates.h"
#include "next_bit_tests.h"

int main(void) {
    printf("=== PRG Example 1: Yao Theorem Verification ===\n\n");

    /* Demonstrate Yao's Theorem: next-bit unpredictability -> pseudorandomness */
    printf("1. Statistical Distance Demonstration:\n");
    double p[] = {0.25, 0.25, 0.25, 0.25};
    double q[] = {0.4, 0.2, 0.1, 0.3};
    double dist = prg_statistical_distance(p, q, 4);
    printf("   delta(p,q) = %.4f\n\n", dist);

    printf("2. Hybrid Argument Bound:\n");
    for (size_t m = 1; m <= 8; m *= 2) {
        double bound = prg_hybrid_argument_bound(m, 0.01);
        printf("   m=%zu bits: adv_dist <= %.4f\n", m, bound);
    }
    printf("\n");

    printf("3. BPP Error Reduction:\n");
    double deltas[] = {0.1, 0.01, 0.001, 1e-6};
    for (int i = 0; i < 4; i++) {
        size_t k = prg_bpp_trials_for_error(deltas[i]);
        printf("   delta=%.0e => trials=%zu\n", deltas[i], k);
    }
    printf("\n");

    printf("4. NW Generator Parameters:\n");
    printf("   d=10, r=2 => m <= %zu bits\n", prg_nw_output_length(10, 2));
    printf("   m=1000, k=8, r=2 => d >= %zu bits\n", prg_nw_design_size(1000, 8, 2));
    printf("\n");

    printf("5. IW97 BPP=P Bound:\n");
    double time_bound = prg_iw97_bound(20, 0.1, 100);
    printf("   m=20, eps=0.1, n=100 => deterministic time ~ %.2e\n", time_bound);
    printf("\n");

    printf("6. PRG Security Analysis:\n");
    double sec = prg_security_level(128, 1e-10);
    printf("   128-bit seed, adv=1e-10 => security = %.1f bits\n", sec);
    size_t min_seed = prg_min_seed_length(128, 1e-8);
    printf("   Target 128-bit security, adv=1e-8 => min seed = %zu bits\n", min_seed);
    printf("\n");

    printf("7. Negligibility Check:\n");
    printf("   adv=2^-50, n=64 => %s\n",
           prg_is_negligible(pow(2.0, -50.0), 64) ? "negligible" : "not negligible");
    printf("   adv=0.1, n=128 => %s\n",
           prg_is_negligible(0.1, 128) ? "negligible" : "not negligible");

    printf("\n=== Example Complete ===\n");
    return 0;
}
