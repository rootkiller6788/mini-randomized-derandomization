/******************************************************************************
 * bench_circuit_lb.c ¡ª Performance Benchmarks for Hardness vs Randomness
 *
 * Benchmarks core operations: circuit evaluation, lower bound computation,
 * truth table generation, derandomization analysis.
 *
 * Measures operations/second for key algorithms.
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include "circuit_lower_bounds.h"
#include "hardness_randomness.h"
#include "derandomization_via_hardness.h"
#include "worst_to_average.h"

static double time_now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

#define BENCH(name, N, code) do { \
    double t0 = time_now(); \
    for (size_t _i = 0; _i < (N); _i++) { code; } \
    double t1 = time_now(); \
    double ops = (double)(N) / (t1 - t0); \
    printf("  %-30s %10.0f ops/s\n", name, ops); \
} while(0)

int main(void) {
    printf("=== Hardness vs Randomness Benchmarks ===\n\n");
    size_t iterations = 100000;
    /* Circuit benchmarks */
    BENCH("circuit_create(8)", iterations,
        { BooleanCircuit *c = clb_circuit_create(8); clb_circuit_free(c); });
    /* Lower bound computation */
    BENCH("shannon_counting(8)", iterations,
        { CircuitLB lb = clb_shannon_counting(8); (void)lb; });
    BENCH("hastad_switching(16)", iterations/10,
        { CircuitLB lb = clb_hastad_switching(16, 3); (void)lb; });
    BENCH("razborov_clique(16)", iterations/10,
        { CircuitLB lb = clb_razborov_clique(16); (void)lb; });
    /* HR computations */
    BENCH("hr_hardness_from_circuit_lb", iterations,
        { double h = hr_hardness_from_circuit_lb(16, 1024); (void)h; });
    BENCH("hr_seed_len(WEAK,100)", iterations,
        { size_t s = hr_seed_len(100, HR_WEAK); (void)s; });
    /* WTA computations */
    BENCH("wta_yao_xor_lemma", iterations,
        { WTAConversion c = wta_yao_xor_lemma(0.1, 16, 5); (void)c; });
    BENCH("wta_sample_complexity", iterations,
        { size_t s = wta_sample_complexity(16, 0.1, 0.05); (void)s; });
    /* Circuit evaluation */
    BooleanCircuit *ec = clb_circuit_create(8);
    for (size_t i = 8; i < 16; i++) {
        size_t g = clb_circuit_add_gate(ec, GATE_AND, i-8, i-7, false);
        if (i == 15) clb_circuit_set_output(ec, g);
    }
    bool x[8] = {true, false, true, false, true, false, true, false};
    bool result[1];
    BENCH("circuit_evaluate(8)", iterations,
        { clb_circuit_evaluate(ec, x, result); });
    clb_circuit_free(ec);
    printf("\n=== Benchmarks Complete ===\n");
    return 0;
}