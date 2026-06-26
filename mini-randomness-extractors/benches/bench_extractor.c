/**
 * benches/bench_extractor.c -- Extractor construction benchmarks
 *
 * Compile: cc -I../include -O2 -o bench_ext bench_extractor.c ../src/extractor_core.c ../src/extractor_constructions.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../include/extractor_core.h"
#include "../include/extractor_constructions.h"

static double get_time_sec(void) {
    return (double)clock() / CLOCKS_PER_SEC;
}

int main(void) {
    printf("Extractor Construction Benchmarks\n");
    printf("==================================\n\n");

    double start, elapsed;
    int iters = 10000;

    /* LHL extractor construction */
    start = get_time_sec();
    for (int i = 0; i < iters; i++) {
        SeededExt se = ext_leftover_hash(8, 4);
        bool src[8] = {0}, out[4] = {0};
        se.ext(src, 0, out, 4);
    }
    elapsed = get_time_sec() - start;
    printf("LHL extractor:       %d in %.4f s (%.1f us)\n",
           iters, elapsed, elapsed/iters*1e6);

    /* Trevisan extractor construction */
    start = get_time_sec();
    for (int i = 0; i < iters/10; i++) {
        SeededExt se = ext_trevisan(8, 4, 0.1);
        bool src[8] = {0}, out[4] = {0};
        se.ext(src, i, out, 4);
    }
    elapsed = get_time_sec() - start;
    printf("Trevisan extractor:  %d in %.4f s (%.1f us)\n",
           iters/10, elapsed, elapsed/(iters/10)*1e6);

    /* CG 2-source extractor */
    start = get_time_sec();
    for (int i = 0; i < iters; i++) {
        TwoSourceExt e = ext_two_source_chor_goldreich(8);
        bool x[8] = {0}, y[8] = {0}, out[4] = {0};
        ext_cg_eval(&e, x, y, out);
    }
    elapsed = get_time_sec() - start;
    printf("Chor-Goldreich 2-src: %d in %.4f s (%.1f us)\n",
           iters, elapsed, elapsed/iters*1e6);

    /* Carter-Wegman hash */
    UniversalHash uh = uh_create_carter_wegman(8, 4);
    iters = 500000;
    start = get_time_sec();
    for (int i = 0; i < iters; i++)
        uh_eval(&uh, i % uh.num_a, (i*17) % uh.num_b, i & 0xFF);
    elapsed = get_time_sec() - start;
    printf("Carter-Wegman hash:  %d in %.4f s (%.1f ns)\n",
           iters, elapsed, elapsed/iters*1e9);
    uh_free(&uh);

    /* Disperser eval */
    Disperser d = disp_create(6, 4, 4);
    bool src[6] = {0}, seed[4] = {0}, out[4] = {0};
    iters = 100000;
    start = get_time_sec();
    for (int i = 0; i < iters; i++)
        disp_eval(&d, src, 6, seed, out);
    elapsed = get_time_sec() - start;
    printf("Disperser eval:      %d in %.4f s (%.1f ns)\n",
           iters, elapsed, elapsed/iters*1e9);

    return 0;
}
