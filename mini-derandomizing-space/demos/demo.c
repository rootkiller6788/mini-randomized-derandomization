/**
 * demo.c - Interactive demo of derandomizing space concepts
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "space_derand.h"
#include "logspace_derand.h"
#include "nisan_prg.h"
#include "reingold_zigzag.h"

int main(void) {
    printf("=== mini-derandomizing-space Interactive Demo ===\n\n");
    printf("Classes: %s, %s, %s, %s, %s\n",
           sd_class_name(SC_L), sd_class_name(SC_NL), sd_class_name(SC_RL),
           sd_class_name(SC_BPL), sd_class_name(SC_PSPACE));
    printf("Known containments: L subset NL subset PSPACE\n");
    printf("Savitch bound(10): %zu\n", sd_savitch_space_bound(10));
    printf("NL = coNL: %s\n",
           sd_class_contains(SC_CO_NL, SC_NL) ? "TRUE (Immerman-Szelepcsenyi 1987)" : "FALSE");
    printf("SL = L: %s\n",
           sd_class_contains(SC_SL, SC_L) ? "TRUE (Reingold 2008)" : "FALSE");
    printf("Nisan seed(5,128): %zu\n", sd_nisan_seed_len(5, 128, 0.1));
    RegularGraph *g = sd_cycle_graph(8);
    if (g) {
        SpectralData sd = sd_compute_spectral(g);
        printf("C8 spectral gap: %.6f\n", sd.spectral_gap);
        sd_regular_graph_free(g);
    }
    printf("\n=== Demo Complete ===\n");
    return 0;
}
