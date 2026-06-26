/**
 * @file nw_core.c
 * @brief Core implementation: design operations, PRG evaluation, utilities.
 *
 * Reference:
 *   Nisan & Wigderson, "Hardness vs Randomness", JCSS 49(2), 1994
 *   Arora & Barak, "Computational Complexity: A Modern Approach", Ch 20
 */

#include "nw_core.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* =============================================
 * L1: NWDesign Construction / Destruction
 * ============================================= */

NWDesign *nw_design_create(size_t k, size_t m, size_t l, size_t intersect_bound) {
    if (m == 0 || l == 0 || l > k) return NULL;
    NWDesign *d = (NWDesign *)malloc(sizeof(NWDesign));
    if (!d) return NULL;
    d->k = k; d->m = m; d->l = l; d->intersect_bound = intersect_bound;
    d->sets = (size_t **)malloc(m * sizeof(size_t *));
    if (!d->sets) { free(d); return NULL; }
    for (size_t i = 0; i < m; i++) {
        d->sets[i] = (size_t *)calloc(l, sizeof(size_t));
        if (!d->sets[i]) {
            for (size_t j = 0; j < i; j++) free(d->sets[j]);
            free(d->sets); free(d); return NULL;
        }
    }
    return d;
}

void nw_design_free(NWDesign *d) {
    if (!d) return;
    if (d->sets) {
        for (size_t i = 0; i < d->m; i++) free(d->sets[i]);
        free(d->sets);
    }
    free(d);
}

static int _cmp_size(const void *a, const void *b) {
    size_t va = *(const size_t *)a, vb = *(const size_t *)b;
    return (va > vb) - (va < vb);
}

static size_t _sorted_intersect(const size_t *a, const size_t *b, size_t len) {
    size_t i = 0, j = 0, count = 0;
    while (i < len && j < len) {
        if (a[i] < b[j]) i++;
        else if (a[i] > b[j]) j++;
        else { count++; i++; j++; }
    }
    return count;
}

bool nw_design_verify(const NWDesign *d) {
    if (!d || !d->sets) return false;
    size_t **sorted = (size_t **)malloc(d->m * sizeof(size_t *));
    if (!sorted) return false;
    for (size_t i = 0; i < d->m; i++) {
        sorted[i] = (size_t *)malloc(d->l * sizeof(size_t));
        if (!sorted[i]) {
            for (size_t j = 0; j < i; j++) free(sorted[j]);
            free(sorted); return false;
        }
        memcpy(sorted[i], d->sets[i], d->l * sizeof(size_t));
        qsort(sorted[i], d->l, sizeof(size_t), _cmp_size);
    }
    bool ok = true;
    for (size_t i = 0; i < d->m && ok; i++)
        for (size_t j = i + 1; j < d->m && ok; j++)
            if (_sorted_intersect(sorted[i], sorted[j], d->l) > d->intersect_bound)
                ok = false;
    for (size_t i = 0; i < d->m; i++) free(sorted[i]);
    free(sorted);
    return ok;
}

size_t nw_design_max_intersection(const NWDesign *d) {
    if (!d || !d->sets || d->m < 2) return 0;
    size_t **sorted = (size_t **)malloc(d->m * sizeof(size_t *));
    for (size_t i = 0; i < d->m; i++) {
        sorted[i] = (size_t *)malloc(d->l * sizeof(size_t));
        memcpy(sorted[i], d->sets[i], d->l * sizeof(size_t));
        qsort(sorted[i], d->l, sizeof(size_t), _cmp_size);
    }
    size_t max_inter = 0;
    for (size_t i = 0; i < d->m; i++)
        for (size_t j = i + 1; j < d->m; j++) {
            size_t inter = _sorted_intersect(sorted[i], sorted[j], d->l);
            if (inter > max_inter) max_inter = inter;
        }
    for (size_t i = 0; i < d->m; i++) free(sorted[i]);
    free(sorted);
    return max_inter;
}

void nw_design_print(const NWDesign *d, FILE *fp) {
    if (!d) { fprintf(fp, "NWDesign(NULL)\n"); return; }
    fprintf(fp, "NWDesign(k=%zu, m=%zu, l=%zu, bound=%zu)\n",
            d->k, d->m, d->l, d->intersect_bound);
    size_t show = (d->m < 10) ? d->m : 5;
    for (size_t i = 0; i < show; i++) {
        fprintf(fp, "  S_%zu = {", i);
        for (size_t j = 0; j < d->l && j < 10; j++)
            fprintf(fp, "%zu%s", d->sets[i][j], (j+1 < d->l) ? ", " : "");
        if (d->l > 10) fprintf(fp, "...");
        fprintf(fp, "}\n");
    }
    if (d->m > 10) fprintf(fp, "  ... (%zu more sets)\n", d->m - show);
}

/* =============================================
 * L2: PRG Construction and Evaluation
 * ============================================= */

NWPRG *nw_prg_construct(NWDesign *design, HardFunction *hf) {
    if (!design || !hf) return NULL;
    if (design->l != hf->n) return NULL;
    NWPRG *prg = (NWPRG *)malloc(sizeof(NWPRG));
    if (!prg) return NULL;
    prg->seed_len = design->k;
    prg->output_len = design->m;
    prg->stretch = (design->k > 0) ? (double)design->m / (double)design->k : 0.0;
    prg->design = design;
    prg->hard_fn = hf;
    return prg;
}

void nw_prg_evaluate(const NWPRG *prg, const bool *seed, bool *output) {
    if (!prg || !seed || !output) return;
    bool *restricted = (bool *)malloc(prg->design->l * sizeof(bool));
    if (!restricted) return;
    for (size_t i = 0; i < prg->output_len; i++) {
        for (size_t j = 0; j < prg->design->l; j++) {
            size_t idx = prg->design->sets[i][j];
            restricted[j] = (idx < prg->seed_len) ? seed[idx] : false;
        }
        output[i] = prg->hard_fn->evaluate(restricted, prg->design->l);
    }
    free(restricted);
}

bool nw_prg_evaluate_bit(const NWPRG *prg, const bool *seed, size_t index) {
    if (!prg || !seed || index >= prg->output_len) return false;
    bool *restricted = (bool *)malloc(prg->design->l * sizeof(bool));
    if (!restricted) return false;
    for (size_t j = 0; j < prg->design->l; j++) {
        size_t idx = prg->design->sets[index][j];
        restricted[j] = (idx < prg->seed_len) ? seed[idx] : false;
    }
    bool result = prg->hard_fn->evaluate(restricted, prg->design->l);
    free(restricted);
    return result;
}

void nw_prg_free(NWPRG *prg) { if (prg) free(prg); }

/* =============================================
 * L1/L2: HardFunction Operations
 * ============================================= */

HardFunction *nw_hard_function_create(
    size_t n, bool (*evaluate)(const bool *, size_t),
    double hardness, double epsilon, const char *desc)
{
    HardFunction *hf = (HardFunction *)malloc(sizeof(HardFunction));
    if (!hf) return NULL;
    hf->n = n; hf->evaluate = evaluate;
    hf->hardness = hardness; hf->epsilon = epsilon;
    hf->description = desc ? strdup(desc) : NULL;
    return hf;
}

void nw_hard_function_free(HardFunction *hf) {
    if (!hf) return;
    free(hf->description);
    free(hf);
}

/* =============================================
 * L4: Yao XOR Lemma
 * ============================================= */

double nw_yao_xor_amplify(double epsilon, size_t k, double delta) {
    if (epsilon >= 0.5) return 0.5;
    if (epsilon <= 0.0) return 0.0;
    double advantage = 1.0 - 2.0 * epsilon;
    double amplified_adv = pow(advantage, (double)k);
    double new_epsilon = (1.0 - amplified_adv) / 2.0 + delta;
    if (new_epsilon > 0.5) new_epsilon = 0.5;
    if (new_epsilon < 0.0) new_epsilon = 0.0;
    return new_epsilon;
}

double nw_amplified_circuit_size(double base_size, double epsilon, size_t k, double delta) {
    if (epsilon <= 0 || delta <= 0) return base_size;
    double factor = delta * delta * epsilon * epsilon / (16.0 * (double)(k * k));
    double amplified = base_size * factor;
    if (amplified < 1.0) amplified = 1.0;
    return amplified;
}

/* =============================================
 * L5: Seed Length / Output Length
 * ============================================= */

size_t nw_seed_length(size_t m, double hardness) {
    if (m == 0 || hardness <= 0) return 0;
    if (m == 1) return 1;
    double log_S = log2(hardness);
    if (log_S <= 0) log_S = 1.0;
    size_t l = (size_t)ceil(log_S);
    if (l < 2) l = 2;
    size_t k_base = l * l;
    size_t log_m = nw_ceil_log2(m);
    size_t k = (k_base > log_m) ? k_base : log_m;
    if (k < l) k = l;
    return k;
}

size_t nw_output_length(size_t k, double hardness) {
    if (k == 0 || hardness <= 0) return 0;
    size_t l = (size_t)ceil(sqrt((double)k));
    if (l < 2) l = 2;
    size_t m = (size_t)1 << (l / 2);
    if (m > (size_t)1 << 30) m = (size_t)1 << 30;
    if (m < 2) m = 2;
    return m;
}

/* =============================================
 * L5: NW Theorem Verification
 * ============================================= */

bool nw_theorem_verify(const NWPRG *prg, const NextBitPredictor *predictor, double confidence) {
    (void)predictor; (void)confidence;
    if (!prg || !prg->design) return false;
    size_t log_m = nw_ceil_log2(prg->output_len);
    if (log_m == 0) log_m = 1;
    return (prg->design->intersect_bound <= log_m);
}

NextBitPredictor *nw_hybrid_argument(const NWPRG *prg, const Distinguisher *distinguisher) {
    (void)prg; (void)distinguisher;
    return NULL;
}

/* =============================================
 * L3: Boolean Array Utilities
 * ============================================= */

void nw_int_to_seed(size_t n, uint64_t value, bool *seed) {
    for (size_t i = 0; i < n && i < 64; i++)
        seed[i] = (value >> i) & 1ULL;
    for (size_t i = 64; i < n; i++)
        seed[i] = false;
}

uint64_t nw_seed_to_int(size_t n, const bool *seed) {
    uint64_t result = 0;
    size_t limit = (n < 64) ? n : 64;
    for (size_t i = 0; i < limit; i++)
        if (seed[i]) result |= (1ULL << i);
    return result;
}

void nw_seed_restrict(size_t k, const bool *seed, const size_t *subset, size_t subset_sz, bool *result) {
    for (size_t i = 0; i < subset_sz; i++) {
        size_t idx = subset[i];
        result[i] = (idx < k) ? seed[idx] : false;
    }
}

void nw_bool_xor(const bool *a, const bool *b, bool *out, size_t len) {
    for (size_t i = 0; i < len; i++)
        out[i] = a[i] ^ b[i];
}

size_t nw_hamming_weight(const bool *arr, size_t len) {
    size_t count = 0;
    for (size_t i = 0; i < len; i++)
        if (arr[i]) count++;
    return count;
}

double nw_statistical_distance(const double *dist1, const double *dist2, size_t n) {
    size_t sz = (size_t)1 << n;
    if (n > 20 || sz == 0 || sz > ((size_t)1 << 20)) return 1.0;
    double total = 0.0;
    for (size_t i = 0; i < sz; i++)
        total += fabs(dist1[i] - dist2[i]);
    return total / 2.0;
}

void nw_random_bits(bool *buf, size_t n) {
    static bool seeded = false;
    if (!seeded) { srand((unsigned int)time(NULL)); seeded = true; }
    for (size_t i = 0; i < n; i++)
        buf[i] = (rand() & 1) ? true : false;
}

size_t nw_ceil_log2(size_t n) {
    if (n <= 1) return 0;
    size_t result = 0, val = n - 1;
    while (val > 0) { val >>= 1; result++; }
    return result;
}
