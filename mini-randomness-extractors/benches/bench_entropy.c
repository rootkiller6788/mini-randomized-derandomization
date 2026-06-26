/**
 * benches/bench_entropy.c -- Performance benchmarks for entropy computation
 *
 * Measures throughput of min-entropy and statistical distance computations.
 * Compile: cc -I../include -O2 -o bench_entropy bench_entropy.c ../src/extractor_core.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../include/extractor_core.h"

static double get_time_sec(void) {
    return (double)clock() / CLOCKS_PER_SEC;
}

int main(void) {
    printf("Entropy Computation Benchmarks\n");
    printf("==============================\n\n");

    size_t nvals = 256; /* 8-bit distribution */
    size_t *freq = (size_t*)malloc(nvals * sizeof(size_t));
    for (size_t i = 0; i < nvals; i++) freq[i] = (i * 7919 + 104729) % 100 + 1;
    size_t total = 0;
    for (size_t i = 0; i < nvals; i++) total += freq[i];

    int iterations = 100000;
    double start, elapsed;

    /* Benchmark: min-entropy */
    start = get_time_sec();
    for (int i = 0; i < iterations; i++)
        ext_min_entropy_from_freq(freq, nvals, total);
    elapsed = get_time_sec() - start;
    printf("Min-entropy:     %d iters in %.4f s (%.1f us/iter)\n",
           iterations, elapsed, elapsed/iterations*1e6);

    /* Benchmark: Shannon entropy */
    start = get_time_sec();
    for (int i = 0; i < iterations; i++)
        ext_shannon_entropy(freq, nvals, total);
    elapsed = get_time_sec() - start;
    printf("Shannon entropy: %d iters in %.4f s (%.1f us/iter)\n",
           iterations, elapsed, elapsed/iterations*1e6);

    /* Benchmark: statistical distance */
    size_t *freq2 = (size_t*)malloc(nvals * sizeof(size_t));
    for (size_t i = 0; i < nvals; i++)
        freq2[i] = (i * 31337 + 528459) % 100 + 1;
    size_t total2 = 0;
    for (size_t i = 0; i < nvals; i++) total2 += freq2[i];

    start = get_time_sec();
    for (int i = 0; i < iterations; i++)
        ext_statistical_distance(freq, freq2, nvals, total, total2);
    elapsed = get_time_sec() - start;
    printf("Stat distance:   %d iters in %.4f s (%.1f us/iter)\n",
           iterations, elapsed, elapsed/iterations*1e6);

    /* Benchmark: extractor evaluation */
    Extractor e = ext_create(8, 4, 4, 0.1);
    bool src[8] = {1,0,1,0,1,0,1,0};
    bool seed[4] = {0,0,1,1};
    bool out[4] = {0};
    iterations = 500000;
    start = get_time_sec();
    for (int i = 0; i < iterations; i++)
        ext_eval(&e, src, seed, out);
    elapsed = get_time_sec() - start;
    printf("Extractor eval:  %d iters in %.4f s (%.1f ns/iter)\n",
           iterations, elapsed, elapsed/iterations*1e9);

    free(freq);
    free(freq2);
    return 0;
}
