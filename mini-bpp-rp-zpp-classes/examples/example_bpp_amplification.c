/******************************************************************************
 * example_bpp_amplification.c — Demonstrate BPP Error Amplification
 *
 * Shows how a BPP machine with base error 1/3 can be amplified to
 * arbitrarily small error using independent repetitions and majority vote.
 *
 * L2: Core Concepts — error probability amplification
 * L5: Algorithms — BPP simulation with Chernoff bound guarantees
 *
 * Reference: Arora-Barak §7.1.3, Lemma 7.9
 ******************************************************************************/

#include "randomized_classes.h"
#include "chernoff_bounds.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* Simulated BPP decider with controllable error rate */
static bool biased_decider(const char *input, size_t len, RandomSource *rs) {
    (void)input; (void)len;
    /* 65% chance of correct answer (error = 0.35, close to 1/3) */
    return random_source_uniform(rs, 100) < 65;
}

int main(void) {
    printf("=== BPP Error Amplification Demo ===\n\n");

    printf("A BPP machine with base error 1/3 is amplified to error 2^{-k}\n");
    printf("by taking majority vote over m = 18k independent runs.\n");
    printf("(Chernoff bound: Pr[majority wrong] <= exp(-m/18) = 2^{-k})\n\n");

    printf("%-4s | %-6s | %-12s | %-16s | %s\n",
           "k", "m=18k", "accepts", "confidence", "decision");
    printf("-----+--------+--------------+------------------+----------\n");

    for (size_t k = 1; k <= 10; k++) {
        RandomizedDecision d = bpp_amplify("test", 4, biased_decider, k, 10000);
        printf("%4zu | %6zu | %6zu/%6zu | %.12f | %s\n",
               k, d.trials,
               d.accept_count, d.trials,
               d.confidence,
               d.accepted ? "ACCEPT" : "REJECT");
    }

    printf("\nTheoretical vs empirical comparison:\n");
    printf("  k=1: Theory 2^{-1}=0.5;  Confidence=%.6f\n",
           1.0 - pow(2.0, -1.0));
    printf("  k=5: Theory 2^{-5}=%.5f; Confidence=%.6f\n",
           pow(2.0, -5.0),
           1.0 - pow(2.0, -5.0));
    printf("  k=10: Theory 2^{-10}=%.7f; Confidence=%.7f\n",
           pow(2.0, -10.0),
           1.0 - pow(2.0, -10.0));

    printf("\nComparison RP vs BPP amplification:\n\n");

    printf("RP (OR of trials):\n");
    for (size_t k = 1; k <= 8; k *= 2) {
        size_t trials = rp_trials_for_error(pow(2.0, -(double)k));
        printf("  Error <= 2^{-%zu}: need %zu trials\n", k, trials);
    }
    printf("\nBPP (majority vote):\n");
    for (size_t k = 1; k <= 8; k *= 2) {
        size_t trials = bpp_trials_for_error(pow(2.0, -(double)k));
        printf("  Error <= 2^{-%zu}: need %zu trials\n", k, trials);
    }

    return 0;
}
