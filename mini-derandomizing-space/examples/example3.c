/**
 * example3.c - Zig-Zag Product and Expander Graphs
 * =================================================
 * Demonstrates the zig-zag product, graph powering, and spectral
 * analysis of regular graphs for Reingold's construction.
 *
 * Knowledge: L3 (RegularGraph, SpectralData), L5 (zig-zag product),
 * L8 (expander analysis, Ramanujan graphs).
 * Refs: Reingold-Vadhan-Wigderson (2002), Hoory-Linial-Wigderson (2006).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "space_derand.h"
#include "reingold_zigzag.h"

int main(void) {
    printf("=== Zig-Zag Product & Expander Graphs ===\n\n");

    /* Basic graph creation */
    printf("--- Basic Graphs ---\n");
    RegularGraph *C8 = sd_cycle_graph(8);
    RegularGraph *K5 = sd_complete_graph(5);
    if (C8) printf("  C8: n=%zu, degree=%zu\n", C8->n, C8->d);
    if (K5) printf("  K5: n=%zu, degree=%zu\n", K5->n, K5->d);

    /* Spectral analysis */
    printf("\n--- Spectral Analysis ---\n");
    if (C8) {
        SpectralData sd = sd_compute_spectral(C8);
        printf("  C8 spectral gap: %.6f\n", sd.spectral_gap);
        printf("  C8 lambda_2: %.6f\n", sd.lambda2);
        printf("  C8 edge expansion: %.6f\n", sd.edge_expansion);
    }
    if (K5) {
        SpectralData sd = sd_compute_spectral(K5);
        printf("  K5 spectral gap: %.6f\n", sd.spectral_gap);
        printf("  K5 is expander (lambda<0.5): %s\n",
               sd_is_expander(K5, 0.5) ? "yes" : "no");
    }

    /* Zig-zag product */
    printf("\n--- Zig-Zag Product C4 (R) C2 ---\n");
    RegularGraph *C4 = sd_cycle_graph(4);
    RegularGraph *C2 = sd_cycle_graph(2);
    if (C4 && C2) {
        printf("  C4: n=%zu, d=%zu\n", C4->n, C4->d);
        printf("  C2: n=%zu, d=%zu\n", C2->n, C2->d);
        /* Note: zigzag requires G.d == H.n. Here G=C4(d=2), H=C2(n=2) */
        RegularGraph *Z = sd_zigzag_product(C4, C2);
        if (Z) {
            printf("  C4 (R) C2: n=%zu, d=%zu\n", Z->n, Z->d);
            printf("  Expected: n = 4*2 = 8, d = 2^2 = 4\n");
            bool conn = sd_is_connected(Z);
            printf("  Connected: %s\n", conn ? "yes" : "no");
            sd_regular_graph_free(Z);
        }
    }

    /* Graph powering */
    printf("\n--- Graph Power ---\n");
    if (C8) {
        RegularGraph *C8sq = sd_graph_power(C8, 2);
        if (C8sq) {
            printf("  C8^2: n=%zu, d=%zu\n", C8sq->n, C8sq->d);
            size_t shortest = sd_shortest_path(C8, 0, 4);
            size_t shortest2 = sd_shortest_path(C8sq, 0, 4);
            printf("  C8: dist(0,4) = %zu\n", shortest);
            printf("  C8^2: dist(0,4) = %zu\n", shortest2);
            sd_regular_graph_free(C8sq);
        }
    }

    /* Tensor product */
    printf("\n--- Tensor Product ---\n");
    if (C4) {
        RegularGraph *T = sd_tensor_product(C4, C4);
        if (T) {
            printf("  C4 x C4: n=%zu, d=%zu\n", T->n, T->d);
            sd_regular_graph_free(T);
        }
    }

    /* Logspace connectivity network */
    printf("\n--- Logspace Connectivity Network ---\n");
    RegularGraph *net = rz_logspace_connectivity_network(8);
    if (net) {
        printf("  Network: n=%zu, d=%zu\n", net->n, net->d);
        printf("  Connected: %s\n", sd_is_connected(net) ? "yes" : "no");
    }

    /* Ramanujan check */
    printf("\n--- Ramanujan Graphs ---\n");
    if (K5) printf("  K5 is Ramanujan: %s\n", rz_is_ramanujan(K5) ? "yes" : "no");
    if (C8) printf("  C8 is Ramanujan: %s\n", rz_is_ramanujan(C8) ? "yes" : "no");

    /* Alon-Boppana bound */
    printf("\n--- Alon-Boppana Bound ---\n");
    for (size_t d = 3; d <= 10; d++) {
        double bound = rz_alon_boppana_bound(d);
        printf("  d=%zu: bound = %.6f\n", d, bound);
    }

    /* Cleanup */
    sd_regular_graph_free(C8);
    sd_regular_graph_free(K5);
    sd_regular_graph_free(C4);
    sd_regular_graph_free(C2);
    if (net) sd_regular_graph_free(net);

    printf("\n=== Example3 Complete ===\n");
    return 0;
}
