/**
 * @file nw_design.c
 * @brief Combinatorial design constructions for NW PRG.
 *
 * Implements:
 *   L3: Finite field arithmetic for GF(q)
 *   L5: Reed-Solomon based design construction (explicit)
 *   L5: Random design via probabilistic method
 *   L8: Expander-based design construction
 *   L4: Probabilistic analysis of random designs
 *
 * Reference:
 *   Nisan & Wigderson, "Hardness vs Randomness", JCSS 49(2), 1994
 *   Arora & Barak, Ch 20.2.1
 *   Jukna, "Boolean Function Complexity", Ch 11
 *   Hoory, Linial, Wigderson, "Expander Graphs and Their Applications", 2006
 */

#include "nw_designs.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Internal qsort comparator */
static int _cmp_size_internal(const void *a, const void *b) {
    size_t va = *(const size_t *)a, vb = *(const size_t *)b;
    return (va > vb) - (va < vb);
}

/* =============================================
 * L3: Finite Field GF(q) for Reed-Solomon Designs
 *
 * We support small prime fields (q prime, q <= 256).
 * For larger fields, use GF(2^m) representation
 * but the principle is identical.
 * ============================================= */

/* Modular arithmetic in Z_q */
static size_t _mod_add(size_t a, size_t b, size_t q) { return (a + b) % q; }
static size_t _mod_mul(size_t a, size_t b, size_t q) { return (a * b) % q; }
/* _mod_sub reserved for future use in field polynomial division */
#ifdef UNUSED_FN
static size_t _mod_sub(size_t a, size_t b, size_t q) { return (a >= b) ? (a - b) : (q + a - b) % q; }
#endif

/* Evaluate polynomial p(x) = sum_{j=0}^{d} coeff[j] * x^j in GF(q) */
static size_t _poly_eval(const size_t *coeff, size_t d, size_t x, size_t q) {
    size_t result = 0;
    size_t x_pow = 1;
    for (size_t j = 0; j <= d; j++) {
        result = _mod_add(result, _mod_mul(coeff[j], x_pow, q), q);
        x_pow = _mod_mul(x_pow, x, q);
    }
    return result;
}

/**
 * L5: Reed-Solomon Design Construction
 *
 * Fields: F_q = {0, 1, ..., q-1} with q prime.
 * Universe: F_q x F_q = {0, ..., q^2 - 1} mapped as (a,b) -> a*q + b.
 *
 * For each polynomial p of degree <= d over F_q,
 * define S_p = { (x, p(x)) : x in F_q }.
 *
 * Properties:
 *   - Number of sets m = q^{d+1}
 *   - Each set size l = q
 *   - Intersection: |S_p ∩ S_{p'}| <= d
 *     (since p(x) = p'(x) at at most d points when p != p')
 *   - Universe size k = q^2
 *
 * This gives a (q^2, q^{d+1}, q, d)-design.
 */
NWDesign *nw_reed_solomon_design(size_t q, size_t degree) {
    if (q < 2 || q > 256) return NULL;

    /* Universe size k = q^2 */
    size_t k = q * q;

    /* Number of polynomials of degree <= d: q^{d+1} */
    size_t m = 1;
    for (size_t i = 0; i <= degree; i++) {
        if (m > ((size_t)-1) / q) return NULL; /* overflow */
        m *= q;
    }
    /* Practical limit */
    if (m > 100000) m = 100000;

    /* Each set has size q */
    size_t l = q;
    size_t intersect_bound = degree;

    NWDesign *design = nw_design_create(k, m, l, intersect_bound);
    if (!design) return NULL;

    /* Enumerate polynomials: coefficients (c_d, ..., c_0) in base q
     * poly_idx = c_0 + c_1 * q + ... + c_d * q^d */
    size_t *coeff = (size_t *)calloc(degree + 1, sizeof(size_t));
    if (!coeff) { nw_design_free(design); return NULL; }

    for (size_t poly_idx = 0; poly_idx < m; poly_idx++) {
        /* Decode polynomial index to coefficients */
        size_t tmp = poly_idx;
        for (size_t j = 0; j <= degree; j++) {
            coeff[j] = tmp % q;
            tmp /= q;
        }

        /* Generate set: S_{poly} = { (x, p(x)) : x in F_q } */
        for (size_t x = 0; x < q; x++) {
            size_t y = _poly_eval(coeff, degree, x, q);
            design->sets[poly_idx][x] = x * q + y;  /* encode (x,y) as element of [q^2] */
        }

        /* Sort each set for efficient intersection checking */
        qsort(design->sets[poly_idx], l, sizeof(size_t), _cmp_size_internal);
    }

    free(coeff);
    return design;
}


NWDesign *nw_reed_solomon_design_param(size_t k, size_t m, size_t l, size_t bound) {
    /* Try to find q and degree d such that:
     *   q = l
     *   q^2 = k
     *   q^{d+1} >= m
     *   d <= bound
     */
    size_t q = l;
    if (q < 2) return NULL;
    if (q * q > k) return NULL; /* Universe too small */

    /* Find minimal d such that q^{d+1} >= m */
    size_t d = 0;
    size_t sets = q;
    while (sets < m && d <= bound) {
        sets *= q;
        d++;
    }

    if (d > bound) {
        /* Try with larger q (over GF extension) */
        for (q = l; q <= 256; q++) {
            if (q * q > k) break;
            d = 0;
            sets = q;
            while (sets < m && d <= bound) {
                sets *= q;
                d++;
            }
            if (d <= bound) break;
        }
    }

    if (d > bound) return NULL;
    return nw_reed_solomon_design(q, d);
}

/* =============================================
 * L5: Random Design via Probabilistic Method
 *
 * Each of the m l-element subsets of [k] is chosen
 * independently and uniformly at random.
 *
 * Theorem: When k = l^2 and l = O(log m), a random
 * design has intersection <= log m with probability >= 1/2.
 * ============================================= */

NWDesign *nw_random_design(size_t k, size_t m, size_t l) {
    if (k < l || m == 0) return NULL;

    size_t bound = nw_ceil_log2(m);
    if (bound == 0) bound = 1;

    NWDesign *design = nw_design_create(k, m, l, bound);
    if (!design) return NULL;

    /* Actually we should pick random l-subsets. For realism, we use
     * a deterministic pseudo-random selection via linear congruential
     * generators. This is not cryptographically secure but models
     * the probabilistic method construction. */
    unsigned int rng_state = 12345;
    for (size_t i = 0; i < m; i++) {
        /* Reservoir sampling: select l distinct elements from [k] */
        bool *chosen = (bool *)calloc(k, sizeof(bool));
        if (!chosen) { nw_design_free(design); return NULL; }
        for (size_t j = 0; j < l; j++) {
            size_t remaining = k - j;
            rng_state = rng_state * 1103515245 + 12345;
            size_t skip = rng_state % remaining;
            size_t pos = 0;
            size_t count = 0;
            while (pos < k) {
                if (!chosen[pos]) {
                    if (count == skip) break;
                    count++;
                }
                pos++;
            }
            if (pos >= k) pos = k - 1;
            while (chosen[pos] && pos < k) pos++;
            if (pos >= k) {
                /* fallback: find any unchosen */
                for (pos = 0; pos < k && chosen[pos]; pos++) {}
            }
            chosen[pos] = true;
            design->sets[i][j] = pos;
        }
        free(chosen);
        qsort(design->sets[i], l, sizeof(size_t), _cmp_size_internal);
    }

    return design;
}

/**
 * L4: Probabilistic Analysis of Random Designs
 *
 * For a random design with parameters (k, m, l):
 *   P[|S_i ∩ S_j| > t] <= C(l, t+1) * (l/k)^{t+1} * (C(k-t-1, l-t-1) / C(k, l))
 *
 * Using union bound over C(m, 2) pairs:
 *   P[failure] <= C(m, 2) * P[pair fails]
 *
 * Since we sort each set, we compute exact via hypergeometric:
 *   Expected intersection = l^2 / k
 */
double nw_random_design_success_prob(size_t k, size_t m, size_t l, size_t bound) {
    if (l > k) return 0.0;
    if (bound >= l) return 1.0;

    /* For random l-subsets of [k], intersection size follows
     * hypergeometric distribution with mean l^2/k.
     *
     * We use Chernoff bound for hypergeometric:
     *   P[|S_i ∩ S_j| >= mean + δ] <= exp(-2δ^2 / l)
     *
     * Let δ = bound + 1 - mean
     */
    double mean = (double)(l * l) / (double)k;
    double excess = (double)(bound + 1) - mean;
    if (excess <= 0) return 0.0; /* Expected intersection already exceeds bound */

    double single_pair_fail = exp(-2.0 * excess * excess / (double)l);
    if (single_pair_fail >= 1.0) single_pair_fail = 1.0 - 1e-10;

    /* Union bound over C(m,2) pairs */
    double num_pairs = (double)m * (double)(m - 1) / 2.0;
    double union_fail = num_pairs * single_pair_fail;

    double success = 1.0 - union_fail;
    if (success < 0.0) success = 0.0;
    if (success > 1.0) success = 1.0;
    return success;
}

size_t nw_design_trials_needed(size_t k, size_t m, size_t l, size_t bound, double confidence) {
    double p = nw_random_design_success_prob(k, m, l, bound);
    if (p <= 0.0) return (size_t)-1;  /* Infeasible */
    if (p >= confidence) return 1;

    /* trials = ceil(log(1-confidence) / log(1-p)) */
    double log_target = log(1.0 - confidence);
    double log_fail = log(1.0 - p);
    if (log_fail >= 0.0) return (size_t)-1; /* p=0 effectively */

    double trials = ceil(log_target / log_fail);
    return (trials < 1) ? 1 : (size_t)trials;
}

/* =============================================
 * L8: Expander Graph Construction (Margulis)
 *
 * Margulis (1973) gave the first explicit expander.
 * Here we implement a simplified version for n = m^2
 * with degree 8.
 *
 * Vertices are pairs (x,y) in Z_m x Z_m.
 * Neighbors of (x,y):
 *   T1: (x, y)
 *   T2: (x, x+y)
 *   T3: (x, x+y+1)
 *   T4: (x+y, y)
 *   T5: (x+y+1, y)
 *   plus inverses (total degree <= 8)
 *
 * This graph is known to have λ <= 5√2 ≈ 7.07 < 8.
 * ============================================= */

Expander *nw_expander_create_margulis(size_t n, size_t d) {
    (void)d; /* degree is fixed at 8 for Margulis */
    if (n < 4) return NULL;

    /* Find m such that m^2 >= n */
    size_t m = (size_t)ceil(sqrt((double)n));
    n = m * m; /* actual number of vertices */

    Expander *e = (Expander *)malloc(sizeof(Expander));
    if (!e) return NULL;

    e->n = n;
    e->d = 8; /* Margulis expander has degree 8 */
    e->lambda = 5.0 * sqrt(2.0) / 8.0; /* Normalized: λ/d ≈ 0.88 */
    e->adj = (size_t *)calloc(n * 8, sizeof(size_t));
    if (!e->adj) { free(e); return NULL; }

    /* Build adjacency */
    for (size_t x = 0; x < m; x++) {
        for (size_t y = 0; y < m; y++) {
            size_t v = x * m + y;

            /* 8 transformations mod m (Margulis expander) */
            size_t shifts[8][2] = {
                {x, y},
                {x, (x+y) % m},
                {x, (x+y+1) % m},
                {(x+y) % m, y},
                {(x+y+1) % m, y},
                {(x-y+m) % m, y},
                {(x-y-1+m) % m, y},
                {(x+y-1+m) % m, y}
            };

            size_t adj_idx = 0;
            for (size_t s = 1; s <= 8 && adj_idx < 8; s++) {
                if (s < 9) {
                    size_t nx = shifts[s][0];
                    size_t ny = shifts[s][1];
                    if (nx < m && ny < m) {
                        e->adj[v * 8 + adj_idx] = nx * m + ny;
                        adj_idx++;
                    }
                }
            }
            /* Fill remaining with self-loops if needed */
            while (adj_idx < 8) {
                e->adj[v * 8 + adj_idx] = v;
                adj_idx++;
            }
        }
    }

    return e;
}

/**
 * L8: Design from Expander Random Walks
 *
 * Each vertex defines a set: the collection of vertices
 * reachable by length-l non-backtracking random walks.
 * The expander mixing property ensures small intersections.
 */
NWDesign *nw_design_from_expander(const Expander *e, size_t l) {
    if (!e || l == 0) return NULL;

    /* Each vertex gives a set: sample l random steps */
    size_t k = e->n;
    size_t m = e->n; /* one set per vertex */
    size_t bound = nw_ceil_log2(m);
    if (bound < 2) bound = 2;

    NWDesign *design = nw_design_create(k, m, l, bound);
    if (!design) return NULL;

    for (size_t v = 0; v < m; v++) {
        size_t curr = v;
        design->sets[v][0] = curr;
        for (size_t step = 1; step < l; step++) {
            /* Deterministic "random" step: hash-based selection */
            size_t neighbor_idx = (curr * 2654435761ULL + step * 1597334677ULL) % e->d;
            curr = e->adj[curr * e->d + neighbor_idx];
            design->sets[v][step] = curr;
        }
        qsort(design->sets[v], l, sizeof(size_t), _cmp_size_internal);
    }

    return design;
}

double nw_expander_alon_boppana(double d) {
    /* Alon-Boppana: For any infinite family of d-regular graphs,
     * liminf λ >= 2√(d-1).
     * Ramanujan graphs achieve λ <= 2√(d-1) exactly. */
    if (d <= 1) return d;
    return 2.0 * sqrt(d - 1.0);
}

void nw_expander_free(Expander *e) {
    if (!e) return;
    free(e->adj);
    free(e);
}

/* =============================================
 * SetSystem Operations
 * ============================================= */

SetSystem *nw_setsystem_create(size_t n, size_t m, size_t l) {
    SetSystem *ss = (SetSystem *)malloc(sizeof(SetSystem));
    if (!ss) return NULL;
    ss->universe = n;
    ss->num_sets = m;
    ss->set_size = l;
    ss->elements = (size_t *)calloc(m * l, sizeof(size_t));
    if (!ss->elements) { free(ss); return NULL; }
    return ss;
}

size_t nw_setsystem_get(const SetSystem *ss, size_t i, size_t j) {
    if (!ss || i >= ss->num_sets || j >= ss->set_size) return 0;
    return ss->elements[i * ss->set_size + j];
}

void nw_setsystem_set(SetSystem *ss, size_t i, size_t j, size_t val) {
    if (!ss || i >= ss->num_sets || j >= ss->set_size) return;
    ss->elements[i * ss->set_size + j] = val;
}

bool nw_setsystem_verify(const SetSystem *ss, size_t bound) {
    if (!ss || ss->num_sets < 2) return true;
    for (size_t i = 0; i < ss->num_sets; i++) {
        for (size_t j = i + 1; j < ss->num_sets; j++) {
            if (nw_setsystem_intersection(ss, i, j) > bound)
                return false;
        }
    }
    return true;
}

size_t nw_setsystem_intersection(const SetSystem *ss, size_t i, size_t j) {
    if (!ss) return 0;
    size_t count = 0;
    for (size_t a = 0; a < ss->set_size; a++) {
        size_t val_a = ss->elements[i * ss->set_size + a];
        for (size_t b = 0; b < ss->set_size; b++) {
            if (ss->elements[j * ss->set_size + b] == val_a) {
                count++;
                break;
            }
        }
    }
    return count;
}

void nw_setsystem_free(SetSystem *ss) {
    if (!ss) return;
    free(ss->elements);
    free(ss);
}

NWDesign *nw_setsystem_to_design(const SetSystem *ss, size_t bound) {
    if (!ss) return NULL;
    NWDesign *d = nw_design_create(ss->universe, ss->num_sets, ss->set_size, bound);
    if (!d) return NULL;
    for (size_t i = 0; i < ss->num_sets; i++)
        for (size_t j = 0; j < ss->set_size; j++)
            d->sets[i][j] = ss->elements[i * ss->set_size + j];
    return d;
}

SetSystem *nw_design_to_setsystem(const NWDesign *d) {
    if (!d) return NULL;
    SetSystem *ss = nw_setsystem_create(d->k, d->m, d->l);
    if (!ss) return NULL;
    for (size_t i = 0; i < d->m; i++)
        for (size_t j = 0; j < d->l; j++)
            ss->elements[i * d->l + j] = d->sets[i][j];
    return ss;
}

/* =============================================
 * Binomial Coefficients
 * ============================================= */

double nw_binomial(size_t n, size_t k) {
    if (k > n) return 0.0;
    if (k == 0 || k == n) return 1.0;
    if (k > n - k) k = n - k; /* symmetry */

    double result = 1.0;
    for (size_t i = 1; i <= k; i++) {
        result *= (double)(n - k + i);
        result /= (double)i;
        if (result > 1e100) return 1e100; /* overflow guard */
    }
    return result;
}

double nw_log_binomial(size_t n, size_t k) {
    if (k > n) return -INFINITY;
    if (k == 0 || k == n) return 0.0;
    if (k > n - k) k = n - k;

    /* Stirling: ln(C(n,k)) ≈ n·H(k/n) - 0.5·ln(2π·k·(n-k)/n) */
    double p = (double)k / (double)n;
    double entropy = -p * log(p) - (1.0 - p) * log(1.0 - p);
    double result = (double)n * entropy;
    if (k > 0 && (n - k) > 0)
        result -= 0.5 * log(2.0 * M_PI * (double)k * (double)(n - k) / (double)n);
    return result;
}
