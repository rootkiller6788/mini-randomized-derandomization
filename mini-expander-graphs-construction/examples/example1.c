/************************************************************
 * example1.c -- Margulis Expander Construction Demo
 *
 * Demonstrates the explicit construction of an 8-regular
 * Margulis-Gabber-Galil expander graph, computes its spectral
 * properties, and verifies expansion.
 *
 * The Margulis construction was the first explicit infinite
 * family of expanders, proving existence without randomness.
 *
 * Build: make example1
 ************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "expander_core.h"
#include "expander_constructions.h"

int main(void) {
    printf("=== Example 1: Margulis-Gabber-Galil Expander ===\n\n");

    /* Build an 8-regular expander on 25 vertices (5x5 grid) */
    printf("Constructing Margulis expander on n=25 vertices...\n");
    ExpanderGraph *g = exp_margulis_gabber_galil(25);
    if (!g) {
        printf("ERROR: Construction failed (n must be a perfect square)\n");
        return 1;
    }

    printf("Graph: n=%zu vertices, d=%zu-regular\n", g->n, g->d);

    /* Verify regularity */
    if (exp_is_regular(g))
        printf("  Regularity: CONFIRMED (%zu-regular)\n", g->d);

    /* Compute spectral properties */
    printf("\nSpectral Analysis:\n");
    double lam2 = exp_second_eigenvalue(g);
    double gap = exp_spectral_gap(g);
    printf("  lambda_2 (second eigenvalue): %.6f\n", lam2);
    printf("  Spectral gap (d - lambda_2):  %.6f\n", gap);

    /* Alon-Boppana bound */
    double ab_bound = exp_alon_boppana_bound(g->d);
    printf("  Alon-Boppana bound (2*sqrt(d-1)): %.6f\n", ab_bound);
    printf("  Ratio lambda_2 / Alon-Boppana:    %.6f\n",
           lam2 / ab_bound);

    /* Cheeger constant bounds */
    double h_low, h_high;
    exp_cheeger_bounds(g, &h_low, &h_high);
    printf("\nCheeger (Edge Expansion) Bounds:\n");
    printf("  Lower bound: %.6f\n", h_low);
    printf("  Upper bound: %.6f\n", h_high);

    /* Mixing time */
    size_t mix_time = exp_random_walk_mix_time(g, 0.01);
    printf("\nMixing time (eps=0.01): %zu steps\n", mix_time);

    /* Diameter */
    size_t diam = exp_diameter(g);
    printf("Diameter bound: %zu\n", diam);

    /* Connectivity */
    printf("Connected: %s\n", exp_is_connected(g) ? "yes" : "no");
    printf("Bipartite: %s\n", exp_is_bipartite(g) ? "yes" : "no");

    /* Number of edges */
    size_t m = exp_num_edges(g);
    printf("Number of edges: %zu (expected %zu)\n", m, g->n * g->d / 2);

    /* Is it Ramanujan? */
    if (exp_is_ramanujan(g))
        printf("\n>>> This graph is RAMANUJAN (optimal spectral expansion) <<<\n");
    else
        printf("\n>>> Not Ramanujan (lambda_2 > 2*sqrt(d-1)) <<<\n");

    /* Demonstrate random walk */
    printf("\nRandom Walk Simulation:\n");
    size_t pos = 0;
    printf("  Start: vertex %zu", pos);
    for (int i = 0; i < 5; i++) {
        pos = exp_random_walk_step(g, pos);
        printf(" -> %zu", pos);
    }
    printf("\n");

    /* Full summary */
    printf("\n--- Full Graph Summary ---\n");
    exp_print_summary(g, stdout);

    exp_free(g);
    printf("\n=== Example 1 Complete ===\n");
    return 0;
}
