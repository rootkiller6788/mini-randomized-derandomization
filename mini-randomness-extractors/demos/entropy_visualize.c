/**
 * demos/entropy_visualize.c -- Entropy measurement visualization
 *
 * Demonstrates entropy computation on distributions of varying bias.
 * Compile: cc -I../include -o entropy_viz entropy_visualize.c ../src/extractor_core.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/extractor_core.h"

int main(void) {
    printf("Entropy vs. Bias Visualization\n");
    printf("================================\n\n");
    printf("Distribution: 2 outcomes (biased coin)\n");
    printf("n=1, varying Pr[X=1] from 0.5 to 1.0\n\n");
    printf("%-10s %-10s %-10s %-10s\n", "Bias", "H_inf", "H_shan", "H2");
    printf("---------- ---------- ---------- ----------\n");

    for (int pct = 50; pct <= 99; pct++) {
        double bias = pct / 100.0;
        size_t ns = 10000;
        size_t freq[2] = {
            (size_t)(ns * (1.0 - bias)),
            (size_t)(ns * bias)
        };
        double h_inf = ext_min_entropy_from_freq(freq, 2, ns);
        double h_shan = ext_shannon_entropy(freq, 2, ns);
        double h2 = ext_renyi_entropy(freq, 2, ns, 2.0);
        printf("%.2f      %.4f      %.4f      %.4f\n",
               bias, h_inf, h_shan, h2);
    }
    return 0;
}
