#ifndef NEXT_BIT_TESTS_H
#define NEXT_BIT_TESTS_H
#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

typedef enum {
    TEST_FREQUENCY_MONOBIT = 0,
    TEST_FREQUENCY_BLOCK   = 1,
    TEST_RUNS              = 2,
    TEST_LONGEST_RUN       = 3,
    TEST_BINARY_MATRIX     = 4,
    TEST_LINEAR_COMPLEXITY = 5,
    TEST_SERIAL            = 6,
    TEST_APPROX_ENTROPY    = 7,
    TEST_CUMULATIVE_SUMS   = 8,
    TEST_RANDOM_EXCURSIONS = 9,
} NISTTestType;

typedef struct {
    NISTTestType type;
    double p_value;
    bool passed;
} NISTTestResult;

double nist_frequency_monobit(const bool *bits, size_t n);
double nist_frequency_block(const bool *bits, size_t n, size_t block_size);
double nist_runs_test(const bool *bits, size_t n);
double nist_longest_run_ones(const bool *bits, size_t n);
double nist_serial_test(const bool *bits, size_t n, size_t m);
double nist_approximate_entropy(const bool *bits, size_t n, size_t m);
double nist_cumulative_sums(const bool *bits, size_t n, bool forward);
double chi_squared_test(const size_t *observed, const double *expected, size_t buckets);
double kolmogorov_smirnov_test(const double *samples, size_t n);
double binary_entropy(const bool *bits, size_t n);
size_t linear_complexity(const bool *bits, size_t n);
double autocorrelation(const bool *bits, size_t n, size_t lag);
double next_bit_test_battery(const bool *bits, size_t n, NISTTestResult *results, size_t max_tests);
#endif
