/* bench_prg.c — Performance benchmark for NW PRG evaluation */
#include <stdio.h>
#include <time.h>
#include "nw_core.h"
#include "nw_prg.h"
int main(void) {
    printf("NW PRG Performance Benchmark
");
    NWPRG *prg = nw_prg_construct_optimal(1024, 10000.0, NW_PRG_STANDARD);
    if (!prg) { printf("Construction failed.
"); return 1; }
    printf("PRG: seed=%zu, output=%zu
", prg->seed_len, prg->output_len);
    bool *seed = calloc(prg->seed_len, sizeof(bool));
    bool *output = calloc(prg->output_len, sizeof(bool));
    size_t iterations = 10000;
    clock_t start = clock();
    for (size_t i = 0; i < iterations; i++)
        nw_prg_evaluate(prg, seed, output);
    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Evaluations: %zu in %.3fs (%.0f eval/s)
", iterations, elapsed, iterations/elapsed);
    free(seed); free(output);
    nw_prg_free(prg);
    return 0;
}
