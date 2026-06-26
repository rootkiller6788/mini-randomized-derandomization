#ifndef WTA_H
#define WTA_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>

/* Result of worst-case to average-case conversion */
typedef struct {
    size_t n;
    double worst_hardness;
    double avg_hardness;
} WTAConversion;

/* Hardcore set from Impagliazzo's lemma */
typedef struct {
    size_t n;
    double density;
    bool *set;
    double hardcore_hardness;
} HardcoreSet;

WTAConversion wta_yao_xor_lemma(double worst, size_t n, size_t k);
WTAConversion wta_impagliazzo_hardcore(double worst, size_t n);
double wta_average_from_worst(double worst, size_t n);
bool wta_conversion_feasible(double worst, size_t n);
size_t wta_sample_complexity(size_t n, double epsilon, double delta);
double wta_local_list_decode_hardness(double worst, size_t n, size_t list_size);
bool wta_goldreich_levin_applies(size_t n, double hardness);
HardcoreSet *wta_construct_hardcore_set(double worst, size_t n, double density);
void wta_free_hardcore_set(HardcoreSet *hs);
double wta_xor_lemma_bound(double worst, size_t k);
size_t wta_optimal_xor_copies(double worst, double target_avg);

#endif /* WTA_H */