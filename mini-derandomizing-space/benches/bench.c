/**
 * bench.c - Performance benchmarks for derandomizing space operations
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "space_derand.h"
#include "logspace_derand.h"
#include "nisan_prg.h"
#include "reingold_zigzag.h"

static double now_sec(void) { return (double)clock() / (double)CLOCKS_PER_SEC; }

int main(void) {
    printf("=== Benchmarks ===\n");
    double t0, t1;

    /* Benchmark: USTCONN on cycle */
    t0 = now_sec();
    size_t edges[100];
    for (size_t i = 0; i < 50; i++) {
        edges[2*i] = i; edges[2*i+1] = (i+1) % 50;
    }
    bool conn;
    for (int r = 0; r < 100; r++)
        sd_reingold_ustconn(50, edges, 50, 0, 25, &conn);
    t1 = now_sec();
    printf("USTCONN(50-cycle, 100 reps): %.6f sec\n", t1 - t0);

    /* Benchmark: complete graph spectral analysis */
    t0 = now_sec();
    RegularGraph *K = sd_complete_graph(20);
    for (int r = 0; r < 50; r++) {
        SpectralData sd = sd_compute_spectral(K);
        (void)sd;
    }
    t1 = now_sec();
    printf("K20 spectral(50 reps): %.6f sec\n", t1 - t0);
    sd_regular_graph_free(K);

    /* Benchmark: hash family */
    t0 = now_sec();
    for (int r = 0; r < 10; r++) {
        PairwiseHashFamily *phf = sd_hash_family_create(8, 8);
        for (size_t i = 0; i < 100; i++)
            sd_hash_eval(phf, i % phf->num_functions, i);
        sd_hash_family_free(phf);
    }
    t1 = now_sec();
    printf("Hash family(10 creates, 1000 evals): %.6f sec\n", t1 - t0);

    /* Benchmark: Nisan PRG */
    t0 = now_sec();
    SpacePRG *prg = nisan_prg_create(5, 64);
    if (prg) {
        bool seed[32] = {0}, out[128] = {0};
        for (int r = 0; r < 500; r++)
            nisan_prg_eval(prg, seed, prg->seed_len, out, prg->output_len);
        nisan_prg_free(prg);
    }
    t1 = now_sec();
    printf("Nisan PRG eval(500 reps): %.6f sec\n", t1 - t0);

    printf("=== Benchmarks Complete ===\n");
    return 0;
}
