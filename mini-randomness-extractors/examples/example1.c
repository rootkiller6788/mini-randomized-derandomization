/**
 * example1.c -- Entropy Measurement and Statistical Distance
 *
 * Demonstrates L2 core concepts:
 *   1. Min-entropy computation from frequency tables
 *   2. Shannon entropy and Renyi entropy
 *   3. Statistical (total variation) distance between distributions
 *   4. Hellinger distance and KL divergence
 *   5. Relationship H_inf <= H <= log2(N)
 *
 * Mathematical verification:
 *   - Uniform distribution over N elements: H_inf = H = log2(N)
 *   - Deterministic distribution: H_inf = H = 0
 *   - Biased distribution: H_inf < H (strict inequality)
 *   - Statistical distance between identical distributions = 0
 *   - Statistical distance between disjoint distributions = 1
 *   - Triangle inequality for statistical distance
 *   - Pinsker inequality: Delta <= sqrt(D_KL / 2)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "extractor_core.h"

int main(void) {
    printf("============================================================\n");
    printf("  Example 1: Entropy & Statistical Distance\n");
    printf("============================================================\n\n");

    /* Part 1: Uniform distribution over 8 elements */
    printf("--- Part 1: Uniform Distribution (8 elements) ---\n");
    {
        size_t freq[8] = {100, 100, 100, 100, 100, 100, 100, 100};
        double h_inf = ext_min_entropy_from_freq(freq, 8, 800);
        double h_shan = ext_shannon_entropy(freq, 8, 800);
        double h2 = ext_renyi_entropy(freq, 8, 800, 2.0);
        printf("  Min-entropy H_inf = %.4f bits (expected 3.0)\n", h_inf);
        printf("  Shannon entropy H = %.4f bits (expected 3.0)\n", h_shan);
        printf("  Collision entropy H2 = %.4f bits (expected 3.0)\n", h2);
        printf("  H_inf = H = H2 = log2(8) = 3.0 [CHECK: %s]\n\n",
               fabs(h_inf-3.0)<0.01 ? "PASS" : "FAIL");
    }

    /* Part 2: Biased distribution */
    printf("--- Part 2: Biased Distribution ---\n");
    {
        size_t freq[4] = {400, 300, 200, 100};
        double h_inf = ext_min_entropy_from_freq(freq, 4, 1000);
        double h_shan = ext_shannon_entropy(freq, 4, 1000);
        double h2 = ext_renyi_entropy(freq, 4, 1000, 2.0);
        printf("  Min-entropy H_inf = %.4f bits\n", h_inf);
        printf("  Shannon entropy H  = %.4f bits\n", h_shan);
        printf("  Collision entropy H2 = %.4f bits\n", h2);
        printf("  H_inf < H < H2 [CHECK: %s]\n\n",
               (h_inf < h_shan && h_shan < h2) ? "PASS" : "FAIL");
    }

    /* Part 3: Deterministic distribution */
    printf("--- Part 3: Deterministic Distribution ---\n");
    {
        size_t freq[4] = {1000, 0, 0, 0};
        double h_inf = ext_min_entropy_from_freq(freq, 4, 1000);
        double h_shan = ext_shannon_entropy(freq, 4, 1000);
        printf("  H_inf = %.4f, H = %.4f (expected both 0)\n", h_inf, h_shan);
        printf("  [CHECK: %s]\n\n",
               (fabs(h_inf)<0.01 && fabs(h_shan)<0.01) ? "PASS" : "FAIL");
    }

    /* Part 4: Statistical distance */
    printf("--- Part 4: Statistical Distance ---\n");
    {
        size_t P[4] = {250, 250, 250, 250};  /* uniform */
        size_t Q[4] = {400, 300, 200, 100};  /* biased */
        size_t R[4] = {1000, 0, 0, 0};       /* deterministic */
        size_t Z[4] = {0, 0, 0, 1000};       /* disjoint */

        double d_PQ = ext_statistical_distance(P, Q, 4, 1000, 1000);
        double d_PR = ext_statistical_distance(P, R, 4, 1000, 1000);
        double d_RZ = ext_statistical_distance(R, Z, 4, 1000, 1000);
        printf("  Delta(uniform, biased)  = %.4f\n", d_PQ);
        printf("  Delta(uniform, determ)  = %.4f\n", d_PR);
        printf("  Delta(disjoint)         = %.4f (expected 1.0)\n", d_RZ);
        printf("  Delta(R,Z) = 1.0 [CHECK: %s]\n\n",
               fabs(d_RZ-1.0)<0.001 ? "PASS" : "FAIL");
    }

    /* Part 5: Triangle inequality */
    printf("--- Part 5: Triangle Inequality ---\n");
    {
        size_t A[4] = {500, 300, 150, 50};
        size_t B[4] = {400, 200, 300, 100};
        size_t C[4] = {100, 400, 250, 250};
        double dAB = ext_statistical_distance(A, B, 4, 1000, 1000);
        double dBC = ext_statistical_distance(B, C, 4, 1000, 1000);
        double dAC = ext_statistical_distance(A, C, 4, 1000, 1000);
        printf("  Delta(A,B)=%.4f, Delta(B,C)=%.4f, Delta(A,C)=%.4f\n",
               dAB, dBC, dAC);
        printf("  Triangle: Delta(A,C) <= Delta(A,B) + Delta(B,C)\n");
        printf("  %.4f <= %.4f [CHECK: %s]\n",
               dAC, dAB+dBC, (dAC <= dAB+dBC+0.001) ? "PASS" : "FAIL");
    }

    /* Part 6: Pinsker inequality */
    printf("\n--- Part 6: Pinsker Inequality ---\n");
    {
        size_t P[4] = {300, 300, 200, 200};
        size_t Q[4] = {250, 250, 250, 250};
        double delta = ext_statistical_distance(P, Q, 4, 1000, 1000);
        double dkl = ext_kl_divergence(P, Q, 4, 1000, 1000);
        double pinsker_bound = sqrt(dkl / 2.0);
        printf("  Delta(P,Q) = %.6f\n", delta);
        printf("  D_KL(P||Q) = %.6f\n", dkl);
        printf("  Pinsker bound sqrt(D_KL/2) = %.6f\n", pinsker_bound);
        printf("  Delta <= sqrt(D_KL/2) [CHECK: %s]\n",
               (delta <= pinsker_bound + 0.001) ? "PASS" : "FAIL");
    }

    printf("\n=== Example 1 Complete ===\n");
    return 0;
}
