/*
 * polynomial.c — Polynomial Operations over Finite Fields
 *
 * Implements univariate and bivariate polynomial arithmetic over
 * finite fields.  These operations form the algorithmic backbone
 * of Sudan's and Guruswami-Sudan's list-decoding algorithms.
 *
 * Key algorithms:
 *   - Polynomial multiplication (convolution)
 *   - Extended Euclidean algorithm (for GCD, division)
 *   - Horner evaluation over GF(q)
 *   - Chien search for root-finding
 *   - Bivariate polynomial construction for Sudan interpolation
 *   - Hasse derivatives for multiplicity constraints
 */

#include "polynomial.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

/* ==================================================================
 *  Univariate Polynomials
 * ================================================================== */

UnivariatePoly *upoly_zero(void)
{
    UnivariatePoly *p = (UnivariatePoly *)calloc(1, sizeof(UnivariatePoly));
    if (!p) return NULL;
    p->coeffs = NULL;
    p->deg = 0;
    p->capacity = 0;
    return p;
}

UnivariatePoly *upoly_constant(const FiniteField *ff, GFElement c)
{
    (void)ff;
    UnivariatePoly *p = (UnivariatePoly *)calloc(1, sizeof(UnivariatePoly));
    if (!p) return NULL;
    p->capacity = 1;
    p->coeffs = (GFElement *)calloc(1, sizeof(GFElement));
    if (!p->coeffs) { free(p); return NULL; }
    p->coeffs[0] = c;
    p->deg = (c != 0) ? 0 : 0;
    return p;
}

UnivariatePoly *upoly_from_coeffs(GFElement *coeffs, size_t deg, size_t capacity)
{
    UnivariatePoly *p = (UnivariatePoly *)calloc(1, sizeof(UnivariatePoly));
    if (!p) { free(coeffs); return NULL; }
    p->coeffs = coeffs;
    p->deg = deg;
    p->capacity = capacity;
    upoly_trim(p);
    return p;
}

UnivariatePoly *upoly_copy(const UnivariatePoly *p)
{
    if (!p) return NULL;
    UnivariatePoly *cp = (UnivariatePoly *)calloc(1, sizeof(UnivariatePoly));
    if (!cp) return NULL;
    cp->deg = p->deg;
    cp->capacity = p->capacity;
    if (p->capacity > 0) {
        cp->coeffs = (GFElement *)calloc(p->capacity, sizeof(GFElement));
        if (!cp->coeffs) { free(cp); return NULL; }
        memcpy(cp->coeffs, p->coeffs, p->capacity * sizeof(GFElement));
    }
    return cp;
}

void upoly_free(UnivariatePoly *p)
{
    if (!p) return;
    free(p->coeffs);
    free(p);
}

void upoly_trim(UnivariatePoly *p)
{
    if (!p || !p->coeffs) return;
    while (p->deg > 0 && p->coeffs[p->deg] == 0) {
        p->deg--;
    }
}

UnivariatePoly *upoly_add(const FiniteField *ff,
                           const UnivariatePoly *a,
                           const UnivariatePoly *b)
{
    assert(ff);
    size_t max_deg = (a->deg > b->deg) ? a->deg : b->deg;
    size_t cap = max_deg + 1;

    GFElement *coeffs = (GFElement *)calloc(cap, sizeof(GFElement));
    if (!coeffs) return NULL;

    for (size_t i = 0; i <= a->deg; i++) coeffs[i] = a->coeffs[i];
    for (size_t i = 0; i <= b->deg; i++)
        coeffs[i] = ff_add(ff, coeffs[i], b->coeffs[i]);

    return upoly_from_coeffs(coeffs, max_deg, cap);
}

UnivariatePoly *upoly_sub(const FiniteField *ff,
                           const UnivariatePoly *a,
                           const UnivariatePoly *b)
{
    assert(ff);
    size_t max_deg = (a->deg > b->deg) ? a->deg : b->deg;
    size_t cap = max_deg + 1;

    GFElement *coeffs = (GFElement *)calloc(cap, sizeof(GFElement));
    if (!coeffs) return NULL;

    for (size_t i = 0; i <= a->deg; i++) coeffs[i] = a->coeffs[i];
    for (size_t i = 0; i <= b->deg; i++)
        coeffs[i] = ff_sub(ff, coeffs[i], b->coeffs[i]);

    return upoly_from_coeffs(coeffs, max_deg, cap);
}

UnivariatePoly *upoly_mul(const FiniteField *ff,
                           const UnivariatePoly *a,
                           const UnivariatePoly *b)
{
    assert(ff);
    if (!a || !b) return NULL;

    size_t res_deg = a->deg + b->deg;
    size_t cap = res_deg + 2;  /* +1 for safety */

    GFElement *coeffs = (GFElement *)calloc(cap, sizeof(GFElement));
    if (!coeffs) return NULL;

    for (size_t i = 0; i <= a->deg; i++) {
        for (size_t j = 0; j <= b->deg; j++) {
            GFElement prod = ff_mul(ff, a->coeffs[i], b->coeffs[j]);
            coeffs[i + j] = ff_add(ff, coeffs[i + j], prod);
        }
    }

    return upoly_from_coeffs(coeffs, res_deg, cap);
}

GFElement upoly_eval(const FiniteField *ff,
                      const UnivariatePoly *p, GFElement x)
{
    assert(ff);
    if (!p || !p->coeffs) return 0;

    /* Horner's method */
    if (p->deg == 0) return p->coeffs[0];

    GFElement result = p->coeffs[p->deg];
    for (size_t i = p->deg; i > 0; i--) {
        result = ff_mul(ff, result, x);
        result = ff_add(ff, result, p->coeffs[i - 1]);
    }

    return result;
}

UnivariatePoly *upoly_derivative(const FiniteField *ff,
                                  const UnivariatePoly *p)
{
    assert(ff);
    if (!p || p->deg == 0) return upoly_constant(ff, 0);

    size_t new_deg = p->deg - 1;
    GFElement *coeffs = (GFElement *)calloc(new_deg + 1, sizeof(GFElement));
    if (!coeffs) return NULL;

    /* d/dx (c_i·x^i) = i·c_i·x^{i-1} mod ff */
    for (size_t i = 1; i <= p->deg; i++) {
        if (p->coeffs[i] != 0) {
            GFElement factor = i % ff->q;
            coeffs[i - 1] = ff_mul(ff, factor, p->coeffs[i]);
        }
    }

    return upoly_from_coeffs(coeffs, new_deg, new_deg + 1);
}

/* Extended Euclidean Algorithm for polynomials over GF(q) */
UnivariatePoly *upoly_extended_gcd(const FiniteField *ff,
                                    const UnivariatePoly *a,
                                    const UnivariatePoly *b,
                                    UnivariatePoly **s_out,
                                    UnivariatePoly **t_out)
{
    assert(ff);
    UnivariatePoly *old_r = upoly_copy(a);
    UnivariatePoly *r = upoly_copy(b);
    UnivariatePoly *old_s = upoly_constant(ff, 1);
    UnivariatePoly *s = upoly_constant(ff, 0);
    UnivariatePoly *old_t = upoly_constant(ff, 0);
    UnivariatePoly *t = upoly_constant(ff, 1);

    while (r->deg > 0 || r->coeffs[0] != 0) {
        UnivariatePoly *q = NULL, *rem = NULL;
        upoly_divmod(ff, old_r, r, &q, &rem);

        /* old_r, r = r, rem */
        UnivariatePoly *tmp_r = old_r; old_r = r; r = rem;

        /* old_s, s = s, old_s - q*s */
        UnivariatePoly *qs = upoly_mul(ff, q, s);
        UnivariatePoly *tmp_s = upoly_sub(ff, old_s, qs);
        upoly_free(qs); upoly_free(old_s); qs = NULL;
        old_s = s; s = tmp_s;

        /* old_t, t = t, old_t - q*t */
        UnivariatePoly *qt = upoly_mul(ff, q, t);
        UnivariatePoly *tmp_t = upoly_sub(ff, old_t, qt);
        upoly_free(qt); upoly_free(old_t); qt = NULL;
        old_t = t; t = tmp_t;

        upoly_free(q);
    }

    /* Make GCD monic */
    if (old_r->deg > 0 || old_r->coeffs[0] != 0) {
        GFElement lead_inv = ff_inv(ff, old_r->coeffs[old_r->deg]);
        for (size_t i = 0; i <= old_r->deg; i++) {
            old_r->coeffs[i] = ff_mul(ff, old_r->coeffs[i], lead_inv);
        }
        for (size_t i = 0; i <= old_s->deg; i++) {
            old_s->coeffs[i] = ff_mul(ff, old_s->coeffs[i], lead_inv);
        }
        for (size_t i = 0; i <= old_t->deg; i++) {
            old_t->coeffs[i] = ff_mul(ff, old_t->coeffs[i], lead_inv);
        }
    }

    if (s_out) *s_out = old_s; else upoly_free(old_s);
    if (t_out) *t_out = old_t; else upoly_free(old_t);
    upoly_free(r); upoly_free(s); upoly_free(t);

    return old_r;
}

/* Polynomial long division: a = q·b + r, deg(r) < deg(b) */
void upoly_divmod(const FiniteField *ff,
                  const UnivariatePoly *a,
                  const UnivariatePoly *b,
                  UnivariatePoly **q_out,
                  UnivariatePoly **r_out)
{
    assert(ff);

    /* Handle degenerate cases */
    if (!b || b->deg == 0) {
        /* b is constant: divide each coefficient */
        if (b && b->coeffs[0] == 0) {
            if (q_out) *q_out = upoly_zero();
            if (r_out) *r_out = upoly_copy(a);
            return;
        }
        /* Divisor is non-zero constant */
        GFElement inv_b0 = ff_inv(ff, b->coeffs[0]);
        size_t cap = a->deg + 1;
        GFElement *qcoeffs = (GFElement *)calloc(cap, sizeof(GFElement));
        if (qcoeffs) {
            for (size_t i = 0; i <= a->deg; i++)
                qcoeffs[i] = ff_mul(ff, a->coeffs[i], inv_b0);
        }
        if (q_out) *q_out = upoly_from_coeffs(qcoeffs, a->deg, cap);
        else free(qcoeffs);
        if (r_out) *r_out = upoly_zero();
        return;
    }

    /* Standard polynomial long division */
    if (a->deg < b->deg) {
        if (q_out) *q_out = upoly_zero();
        if (r_out) *r_out = upoly_copy(a);
        return;
    }

    size_t q_deg = a->deg - b->deg;
    GFElement *q_coeffs = (GFElement *)calloc(q_deg + 1, sizeof(GFElement));
    GFElement *r_coeffs = (GFElement *)calloc(a->deg + 1, sizeof(GFElement));

    /* Copy a into remainder */
    for (size_t i = 0; i <= a->deg; i++) r_coeffs[i] = a->coeffs[i];

    GFElement inv_lead_b = ff_inv(ff, b->coeffs[b->deg]);

    for (size_t k = a->deg; k >= b->deg; k--) {
        /* Not enough for proper sub — check degree */
        if (k - b->deg > q_deg) break;

        GFElement factor = ff_mul(ff, r_coeffs[k], inv_lead_b);
        q_coeffs[k - b->deg] = factor;

        /* Subtract factor · b · x^{k-deg(b)} from remainder */
        for (size_t j = 0; j <= b->deg; j++) {
            GFElement sub = ff_mul(ff, factor, b->coeffs[j]);
            size_t pos = j + (k - b->deg);
            if (pos <= a->deg) {
                r_coeffs[pos] = ff_sub(ff, r_coeffs[pos], sub);
            }
        }
    }

    UnivariatePoly *q = upoly_from_coeffs(q_coeffs, q_deg, q_deg + 1);
    UnivariatePoly *r = upoly_from_coeffs(r_coeffs, a->deg, a->deg + 1);
    upoly_trim(r);

    if (q_out) *q_out = q; else upoly_free(q);
    if (r_out) *r_out = r; else upoly_free(r);
}

/* Chien search: find all roots of p(x) in GF(q) */
size_t upoly_find_roots(const FiniteField *ff,
                         const UnivariatePoly *p,
                         GFElement **roots_out)
{
    assert(ff);
    if (!p || p->deg == 0) {
        *roots_out = NULL;
        return 0;
    }

    /* Upper bound: at most deg roots */
    GFElement *roots = (GFElement *)malloc(p->deg * sizeof(GFElement));
    if (!roots) { *roots_out = NULL; return 0; }

    size_t count = 0;

    /* Try all elements of GF(q) */
    for (GFElement x = 0; x < ff->q && count < p->deg; x++) {
        GFElement val = upoly_eval(ff, p, x);
        if (ff_eq(ff, val, 0)) {
            roots[count++] = x;
        }
    }

    *roots_out = roots;
    return count;
}

/* ==================================================================
 *  Bivariate Polynomials (Sudan / Guruswami-Sudan core)
 * ================================================================== */

BivariatePoly *bpoly_create(size_t deg_x, size_t deg_y)
{
    BivariatePoly *bp = (BivariatePoly *)calloc(1, sizeof(BivariatePoly));
    if (!bp) return NULL;
    bp->deg_x = deg_x;
    bp->deg_y = deg_y;
    size_t n = (deg_x + 1) * (deg_y + 1);
    bp->coeffs = (GFElement *)calloc(n, sizeof(GFElement));
    if (!bp->coeffs) { free(bp); return NULL; }
    return bp;
}

void bpoly_free(BivariatePoly *bp)
{
    if (!bp) return;
    free(bp->coeffs);
    free(bp);
}

static inline size_t bpoly_index(const BivariatePoly *bp, size_t i, size_t j)
{
    return i * (bp->deg_y + 1) + j;
}

GFElement bpoly_coeff(const BivariatePoly *bp, size_t i, size_t j)
{
    if (!bp || i > bp->deg_x || j > bp->deg_y) return 0;
    return bp->coeffs[bpoly_index(bp, i, j)];
}

void bpoly_set_coeff(BivariatePoly *bp, size_t i, size_t j, GFElement val)
{
    if (!bp || i > bp->deg_x || j > bp->deg_y) return;
    bp->coeffs[bpoly_index(bp, i, j)] = val;
}

GFElement bpoly_eval(const FiniteField *ff,
                      const BivariatePoly *bp,
                      GFElement x0, GFElement y0)
{
    assert(ff);
    if (!bp) return 0;

    GFElement result = 0;
    GFElement x_pow = 1;

    for (size_t i = 0; i <= bp->deg_x; i++) {
        GFElement y_pow = 1;
        for (size_t j = 0; j <= bp->deg_y; j++) {
            GFElement term = bpoly_coeff(bp, i, j);
            if (!ff_eq(ff, term, 0)) {
                GFElement tmp = ff_mul(ff, term, x_pow);
                tmp = ff_mul(ff, tmp, y_pow);
                result = ff_add(ff, result, tmp);
            }
            y_pow = ff_mul(ff, y_pow, y0);
        }
        x_pow = ff_mul(ff, x_pow, x0);
    }

    return result;
}

/* Binomial coefficient modulo field characteristic */
static GFElement binom_mod(const FiniteField *ff, size_t n, size_t k)
{
    if (k > n) return 0;
    if (k == 0 || k == n) return 1;
    /* Use multiplicative formula: C(n,k) = prod_{i=1}^k (n-k+i)/i */
    GFElement result = 1;
    for (size_t i = 1; i <= k; i++) {
        result = ff_mul(ff, result, (n - k + i) % ff->q);
        result = ff_div(ff, result, i % ff->q);
    }
    return result;
}

/* Hasse derivative: D_x^a D_y^b Q evaluated at (x0, y0).
 * Uses the definition: coefficient of x^{i}y^{j} is multiplied
 * by C(i,a)·C(j,b) when computing the derivative. */
GFElement bpoly_hasse_derivative(const FiniteField *ff,
                                  const BivariatePoly *bp,
                                  size_t a, size_t b,
                                  GFElement x0, GFElement y0)
{
    assert(ff);
    if (!bp) return 0;

    GFElement result = 0;

    for (size_t i = a; i <= bp->deg_x; i++) {
        GFElement binom_i = binom_mod(ff, i, a);
        if (ff_eq(ff, binom_i, 0)) continue;

        GFElement x_term = ff_mul(ff, binom_i, ff_pow(ff, x0, i - a));

        for (size_t j = b; j <= bp->deg_y; j++) {
            GFElement binom_j = binom_mod(ff, j, b);
            if (ff_eq(ff, binom_j, 0)) continue;

            GFElement coeff = bpoly_coeff(bp, i, j);
            if (ff_eq(ff, coeff, 0)) continue;

            GFElement y_term = ff_mul(ff, binom_j, ff_pow(ff, y0, j - b));
            GFElement term = ff_mul(ff, coeff, ff_mul(ff, x_term, y_term));
            result = ff_add(ff, result, term);
        }
    }

    return result;
}

/* Factor Q(x,y) into linear factors of the form y - p(x).
 * Implementation: for each candidate root y0, check if Q(x, y0) = 0
 * identically (as a polynomial in x).  If so, find the polynomial
 * p(x) = y0 (constant), or more generally, use linear algebra. */
size_t bpoly_factor_linear_y(const FiniteField *ff,
                              const BivariatePoly *Q,
                              UnivariatePoly ***factors_out)
{
    assert(ff);
    if (!Q) { *factors_out = NULL; return 0; }

    /* For a proper implementation, use the Roth-Ruckenstein root-finding
     * or the Gao-Shokrollahi approach.  Here we implement the basic
     * approach: find all polynomials p(x) of degree < k such that
     * Q(x, p(x)) ≡ 0 by solving linear equations. */

    /* Allocate space for at most deg_y factors (Bezout bound) */
    size_t max_factors = Q->deg_y + 1;
    UnivariatePoly **factors = (UnivariatePoly **)calloc(max_factors, sizeof(UnivariatePoly *));
    if (!factors) { *factors_out = NULL; return 0; }

    size_t count = 0;

    /* For each y ∈ GF(q) (small q), check if y - y0 divides Q */
    /* This is a simplified factorization for demonstration.
     * Production code would use the full Sudan root-finding algorithm. */
    for (GFElement y0 = 0; y0 < ff->q && count < max_factors; y0++) {
        /* Compute Q(x, y0): this gets us a univariate polynomial in x */
        GFElement *qx_coeffs = (GFElement *)calloc(Q->deg_x + 1, sizeof(GFElement));
        if (!qx_coeffs) continue;

        bool all_zero = true;
        for (size_t i = 0; i <= Q->deg_x; i++) {
            GFElement sum = 0;
            GFElement y0_pow_j = 1;
            for (size_t j = 0; j <= Q->deg_y; j++) {
                GFElement c = bpoly_coeff(Q, i, j);
                if (!ff_eq(ff, c, 0)) {
                    sum = ff_add(ff, sum, ff_mul(ff, c, y0_pow_j));
                }
                y0_pow_j = ff_mul(ff, y0_pow_j, y0);
            }
            qx_coeffs[i] = sum;
            if (!ff_eq(ff, sum, 0)) all_zero = false;
        }

        if (all_zero) {
            /* y - y0 is a factor → p(x) = y0 (constant polynomial) */
            UnivariatePoly *candidate = upoly_constant(ff, y0);
            if (candidate) {
                factors[count++] = candidate;
            }
        }

        free(qx_coeffs);
    }

    *factors_out = factors;
    return count;
}

/* ==================================================================
 *  Weighted Degree (for Sudan algorithm parameter selection)
 * ================================================================== */

size_t bpoly_weighted_deg(size_t i, size_t j, size_t k)
{
    /* (1, k-1)-weighted degree: wdeg(x^i·y^j) = i + j·(k-1) */
    return i + j * (k - 1);
}

size_t bpoly_total_weighted_deg(const BivariatePoly *bp, size_t k)
{
    if (!bp) return 0;
    size_t max_wdeg = 0;
    for (size_t i = 0; i <= bp->deg_x; i++) {
        for (size_t j = 0; j <= bp->deg_y; j++) {
            if (!ff_eq(NULL, bpoly_coeff(bp, i, j), 0)) {
                /* Note: ff_eq with NULL is a simplification; in practice
                 * we'd need the field.  But for GF(2) the check works. */
                size_t wdeg = bpoly_weighted_deg(i, j, k);
                if (wdeg > max_wdeg) max_wdeg = wdeg;
            }
        }
    }
    return max_wdeg;
}
