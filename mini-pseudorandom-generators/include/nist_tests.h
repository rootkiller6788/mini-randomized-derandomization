#ifndef NIST_TESTS_H
#define NIST_TESTS_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>

typedef enum { NIST_MONOBIT=0, NIST_BLOCK, NIST_RUNS, NIST_LONGRUN, NIST_RANK, NIST_DFT, NIST_TEMPLATE, NIST_UNIVERSAL, NIST_LINEAR, NIST_SERIAL, NIST_ENTROPY, NIST_CUMSUM, NIST_EXCURSIONS, NIST_EXCURSIONS_VAR } NISTTest;
typedef struct { NISTTest type; double p_value; bool passed; } NISTResult;
double nist_test_monobit(const bool *bits, size_t n, size_t param);
double nist_test_block_freq(const bool *bits, size_t n, size_t param);
double nist_test_runs(const bool *bits, size_t n, size_t param);
double nist_test_longest_run(const bool *bits, size_t n, size_t param);
double nist_erfc(double x);
double nist_normal_cdf(double x);
bool nist_test_suite(const bool *bits, size_t n, NISTResult *results, size_t max_tests);
double binary_entropy_metric(const bool *bits, size_t n);
double autocorrelation_metric(const bool *bits, size_t n, size_t lag);
size_t linear_complexity_metric(const bool *bits, size_t n);

#endif /* NIST_TESTS_H */
