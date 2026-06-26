#ifndef HR_H
#define HR_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include "circuit_lower_bounds.h"

/* Hard function: boolean function with quantified circuit hardness */
typedef struct {
    size_t n;
    double hardness;            /* log of minimum circuit size */
    bool (*fn)(const bool*, size_t);
} HardFn;

/* Pseudorandom generator: stretches k random bits to m bits */
typedef struct {
    size_t seed_len;
    size_t output_len;
    bool (*gen)(const bool*, bool*, size_t, size_t);
} PRG;

/* HR framework parameters */
typedef struct {
    HRLevel hardness_level;
    size_t input_size;
    double epsilon;
    double delta;
} HRParams;

bool hr_implies_prg(const HardFn *hf, HRLevel level, PRG *result);
size_t hr_seed_len(size_t n, HRLevel level);
size_t hr_output_len(size_t n, HRLevel level);
double hr_hardness_from_circuit_lb(size_t n, size_t S);
bool hr_bpp_equals_p_check(HRLevel level, size_t n);
double hr_impagliazzo_wigderson_threshold(size_t n);
bool hr_hardness_amplification_possible(double base, size_t n);
size_t hr_nw_prg_trials(size_t n, double hardness, double epsilon);
const char *hr_level_name(HRLevel level);
bool hr_verify_hardness(const HardFn *hf, size_t samples);
double hr_compute_prg_stretch(size_t seed_len, size_t output_len);
size_t hr_optimal_seed_length(double hardness, double epsilon);
bool hr_construct_nw_design(size_t n, size_t m, size_t **sets_out, size_t *set_sizes_out);

#endif /* HR_H */