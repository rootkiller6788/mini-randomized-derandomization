#ifndef DERAND_H
#define DERAND_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include "circuit_lower_bounds.h"

/* Certificate that a function is hard against circuits */
typedef struct {
    size_t n;
    size_t circuit_sz;
    double hardness;
} HardnessCert;

/* BPP algorithm descriptor */
typedef struct {
    size_t n;
    bool (*alg)(const char*, size_t, const bool*, size_t);
} BPPAlg;

/* Derandomization result */
typedef struct {
    bool result;
    size_t time_steps;
    size_t random_bits_used;
    double confidence;
} DerandResult;

bool derand_bpp_with_hardness(const BPPAlg *alg, const HardnessCert *cert, const char *input, size_t ilen, bool *result);
size_t derand_time_complexity(size_t n, double hardness);
bool derand_verify_certificate(const HardnessCert *cert, size_t n);
HardnessCert derand_cert_from_circuit_lb(size_t n, size_t lower_bound);
double derand_error_prob(size_t n, size_t trials, double hardness);
bool derand_conditional_expectations(const BPPAlg *alg, size_t n, bool *result);
size_t derand_advice_size(size_t n, double hardness);
bool derand_simulate_bpp(const BPPAlg *alg, const char *input, size_t ilen, size_t rlen, DerandResult *out);
bool derand_adleman_derandomize(const BPPAlg *alg, size_t n, bool *advice, size_t advice_len);
double derand_confidence_bound(size_t n, size_t trials);

#endif /* DERAND_H */