/************************************************************
 * example2.c -- Derandomized Error Reduction Demo
 *
 * Demonstrates using expander random walks for error reduction
 * in randomized algorithms. This is a key application of
 * expander graphs: they allow amplifying the success probability
 * of RP/BPP algorithms using exponentially fewer random bits
 * than naive independent repetition.
 *
 * We simulate a Monte Carlo primality test and show how an
 * expander walk reduces error probability.
 *
 * Build: make example2
 ************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "expander_core.h"
#include "expander_constructions.h"
#include "expander_applications.h"

/* Simulated RP algorithm: returns true with probability >= 2/3
 * when input is "yes", false always when input is "no".
 * For demo: simulates a test that errs 1/3 of the time. */
static int call_count = 0;

static bool mock_rp_algorithm(void) {
    call_count++;
    /* Simulate: 2/3 probability of correct answer */
    return (rand() % 3) != 0;
}

int main(void) {
    printf("=== Example 2: Derandomized Error Reduction ===\n\n");

    /* Build a good expander (K_10 is 9-regular, great spectral gap) */
    printf("Building expander graph (K_10)...\n");
    ExpanderGraph *g = exp_create(10, 9);
    for (size_t i = 0; i < 10; i++)
        for (size_t j = i + 1; j < 10; j++)
            exp_add_edge(g, i, j);

    double lam2 = exp_second_eigenvalue(g);
    double lam_norm = lam2 / (double)g->d;
    printf("Graph: n=10, d=9, lambda_norm=%.6f\n", lam_norm);
    printf("Spectral gap: %.6f\n", 1.0 - lam_norm);

    /* Error reduction analysis */
    printf("\n--- Error Reduction Analysis ---\n");
    double base_err = 0.333;  /* 1/3 base error (RP) */
    printf("Base error probability: %.4f\n", base_err);

    for (size_t trials = 10; trials <= 100; trials += 30) {
        double amp_err;
        exp_error_reduction(g, trials, base_err, &amp_err);
        printf("  After %3zu trials: error <= %.10f\n", trials, amp_err);
    }

    /* Run RP amplification */
    printf("\n--- RP Algorithm Amplification ---\n");
    call_count = 0;
    double error;
    bool result = exp_amplify_rp_algorithm(g, mock_rp_algorithm, 50, &error);
    printf("RP calls made: %d\n", call_count);
    printf("Result: %s\n", result ? "ACCEPT" : "REJECT");
    printf("Error bound: %.6f\n", error);

    /* Comparison: bits saved */
    printf("\n--- Random Bits Comparison ---\n");
    size_t k = 100, n = 1024;
    size_t d = 16;
    double saved = exp_random_bits_saved(k, n, d);
    printf("Independent sampling: %.0f bits\n",
           (double)k * log2((double)n));
    printf("Expander walk:       %.0f bits\n",
           log2((double)n) + (double)k * log2((double)d));
    printf("Bits saved:           %.0f bits (%.1f%%)\n",
           saved, 100.0 * saved / ((double)k * log2((double)n)));

    /* Chernoff bound comparison */
    printf("\n--- Chernoff Bound Comparison ---\n");
    double chernoff_exp = exp_chernoff_bound_on_expander(g, 100, 0.5, 0.1);
    double chernoff_ind = exp_chernoff_bound_independent(100, 0.5, 0.1);
    printf("Expander Chernoff bound: %.10f\n", chernoff_exp);
    printf("Independent Chernoff:    %.10f\n", chernoff_ind);
    printf("Ratio (expander/indep):  %.4f\n",
           chernoff_exp / (chernoff_ind + 1e-15));

    printf("\n=== Example 2 Complete ===\n");
    exp_free(g);
    return 0;
}
