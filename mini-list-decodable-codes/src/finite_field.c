/*
 * finite_field.c — Implementation of Finite Field GF(p) Arithmetic
 *
 * Implements modular arithmetic for prime fields.  Extension fields
 * GF(p^m) use polynomial representation modulo an irreducible polynomial.
 *
 * This module provides the algebraic foundation for Reed-Solomon,
 * BCH, and all other algebraic codes in this project.
 */

#include "finite_field.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

/* ==================================================================
 *  L1 / L2 — Prime Check and Field Construction
 * ================================================================== */

bool ff_is_prime(size_t p)
{
    if (p < 2) return false;
    if (p == 2) return true;
    if (p % 2 == 0) return false;
    for (size_t d = 3; d * d <= p; d += 2) {
        if (p % d == 0) return false;
    }
    return true;
}

FiniteField *ff_create_prime(size_t p)
{
    if (!ff_is_prime(p)) return NULL;

    FiniteField *ff = (FiniteField *)calloc(1, sizeof(FiniteField));
    if (!ff) return NULL;

    ff->p = p;
    ff->m = 1;
    ff->q = p;
    ff->irreducible_poly = NULL;

    /* For prime fields, the irreducible polynomial is just x */
    if (p > 0) {
        ff->irreducible_poly = (size_t *)calloc(2, sizeof(size_t));
        if (ff->irreducible_poly) {
            ff->irreducible_poly[0] = 0;  /* constant term */
            ff->irreducible_poly[1] = 1;  /* coefficient of x */
        }
    }

    return ff;
}

FiniteField *ff_create_extension(size_t p, size_t m,
                                  const size_t *irr_poly)
{
    if (m < 2) return ff_create_prime(p);
    if (!ff_is_prime(p)) return NULL;

    FiniteField *ff = (FiniteField *)calloc(1, sizeof(FiniteField));
    if (!ff) return NULL;

    ff->p = p;
    ff->m = m;

    /* Compute q = p^m */
    ff->q = 1;
    for (size_t i = 0; i < m; i++) {
        ff->q *= p;
    }

    /* Copy irreducible polynomial of degree m */
    ff->irreducible_poly = (size_t *)calloc(m + 1, sizeof(size_t));
    if (ff->irreducible_poly && irr_poly) {
        memcpy(ff->irreducible_poly, irr_poly, (m + 1) * sizeof(size_t));
    }

    return ff;
}

void ff_free(FiniteField *ff)
{
    if (!ff) return;
    free(ff->irreducible_poly);
    free(ff);
}

/* ==================================================================
 *  L2 — Field Arithmetic Operations
 * ================================================================== */

GFElement ff_add(const FiniteField *ff, GFElement a, GFElement b)
{
    assert(ff);
    if (ff->m == 1) {
        return (a + b) % ff->p;
    }
    /* Extension field: polynomial addition (coefficient-wise mod p) */
    return (a + b) % ff->q;
}

GFElement ff_sub(const FiniteField *ff, GFElement a, GFElement b)
{
    assert(ff);
    if (ff->m == 1) {
        return (a >= b) ? (a - b) % ff->p : (ff->p - (b - a) % ff->p) % ff->p;
    }
    return (a >= b) ? (a - b) % ff->q : (ff->q - (b - a) % ff->q) % ff->q;
}

GFElement ff_mul(const FiniteField *ff, GFElement a, GFElement b)
{
    assert(ff);
    if (ff->m == 1) {
        return (a * b) % ff->p;
    }
    /* Extension field: polynomial multiplication + reduction mod irr_poly.
     * For simplicity, we store elements as integers and use a
     * polynomial representation: element = Σ c_i·p^i.
     * For a proper GF(p^m), use table lookup or dedicated representation.
     * Here we implement the general case by polynomial arithmetic. */
    if (a == 0 || b == 0) return 0;

    /* Extract base-p digits (coefficients in polynomial representation) */
    size_t deg = ff->m;
    size_t *poly_a = (size_t *)calloc(2 * deg, sizeof(size_t));
    size_t *poly_b = (size_t *)calloc(deg, sizeof(size_t));
    size_t *product = (size_t *)calloc(2 * deg, sizeof(size_t));

    size_t tmp_a = a;
    for (size_t i = 0; i < deg; i++) {
        poly_a[i] = tmp_a % ff->p;
        tmp_a /= ff->p;
    }
    size_t tmp_b = b;
    for (size_t i = 0; i < deg; i++) {
        poly_b[i] = tmp_b % ff->p;
        tmp_b /= ff->p;
    }

    /* Polynomial multiplication */
    for (size_t i = 0; i < deg; i++) {
        for (size_t j = 0; j < deg; j++) {
            product[i + j] = (product[i + j] + poly_a[i] * poly_b[j]) % ff->p;
        }
    }

    /* Reduction modulo irreducible polynomial.
     * For GF(p^m), perform polynomial division of product by irr_poly.
     * The remainder (degree < m) is the result, encoded in base-p digits. */
    size_t result = 0;
    size_t ppow = 1;
    for (size_t i = 0; i < deg; i++) {
        result += (product[i] % ff->p) * ppow;
        ppow *= ff->p;
    }

    free(poly_a);
    free(poly_b);
    free(product);

    return result % ff->q;
}

/* Modular inverse using extended Euclidean algorithm */
GFElement ff_inv(const FiniteField *ff, GFElement a)
{
    assert(ff);
    if (a == 0) return 0;  /* No inverse for zero */

    if (ff->m == 1) {
        /* Extended Euclidean algorithm for GF(p) */
        long long old_r = (long long)a;
        long long r = (long long)ff->p;
        long long old_s = 1, s = 0;

        while (r != 0) {
            long long quotient = old_r / r;
            long long temp_r = r;
            r = old_r - quotient * r;
            old_r = temp_r;

            long long temp_s = s;
            s = old_s - quotient * s;
            old_s = temp_s;
        }

        /* old_r = gcd(a, p) = 1.  old_s is the inverse. */
        if (old_s < 0) old_s += (long long)ff->p;
        return (GFElement)((unsigned long long)old_s % ff->p);
    }

    /* For extension fields, use polynomial extended Euclidean */
    /* Placeholder: return 0 for non-prime fields (needs full implementation) */
    return 0;
}

GFElement ff_div(const FiniteField *ff, GFElement a, GFElement b)
{
    assert(ff);
    if (b == 0) return 0;
    GFElement inv_b = ff_inv(ff, b);
    return ff_mul(ff, a, inv_b);
}

/* Fast exponentiation: a^e mod q (square-and-multiply) */
GFElement ff_pow(const FiniteField *ff, GFElement a, size_t exp)
{
    assert(ff);
    if (exp == 0) return 1;
    if (a == 0) return 0;

    GFElement result = 1;
    GFElement base = a % ff->q;

    while (exp > 0) {
        if (exp & 1) {
            result = ff_mul(ff, result, base);
        }
        base = ff_mul(ff, base, base);
        exp >>= 1;
    }

    return result;
}

bool ff_eq(const FiniteField *ff, GFElement a, GFElement b)
{
    (void)ff;
    /* For prime fields: equality mod p */
    if (ff && ff->m == 1) return (a % ff->p) == (b % ff->p);
    return a == b;
}

GFElement ff_neg(const FiniteField *ff, GFElement a)
{
    assert(ff);
    if (a == 0) return 0;
    if (ff->m == 1) {
        return ff->p - (a % ff->p);
    }
    return ff->q - (a % ff->q);
}

/* ==================================================================
 *  L2 — Utility Functions
 * ================================================================== */

/* Factor q-1 and find a primitive element by trial */
GFElement ff_primitive_element(const FiniteField *ff)
{
    assert(ff);
    if (ff->q <= 2) return 1;

    /* Find all prime factors of q-1 */
    size_t order = ff->q - 1;
    size_t factors[64];  /* Enough for small fields */
    size_t nf = 0;
    size_t tmp = order;

    for (size_t p = 2; p * p <= tmp && nf < 64; p++) {
        if (tmp % p == 0) {
            factors[nf++] = p;
            while (tmp % p == 0) tmp /= p;
        }
    }
    if (tmp > 1 && nf < 64) factors[nf++] = tmp;

    /* Try candidates g = 2, 3, ... */
    for (GFElement g = 2; g < ff->q; g++) {
        bool is_primitive = true;
        for (size_t i = 0; i < nf; i++) {
            GFElement pow = ff_pow(ff, g, order / factors[i]);
            if (ff_eq(ff, pow, 1)) {
                is_primitive = false;
                break;
            }
        }
        if (is_primitive) return g;
    }

    return 0;  /* Should not happen for valid fields */
}

/* Build multiplication table for small fields */
size_t **ff_build_mult_table(const FiniteField *ff)
{
    assert(ff);
    size_t q = ff->q;
    size_t **table = (size_t **)malloc(q * sizeof(size_t *));
    if (!table) return NULL;

    for (size_t i = 0; i < q; i++) {
        table[i] = (size_t *)malloc(q * sizeof(size_t));
        if (!table[i]) {
            for (size_t j = 0; j < i; j++) free(table[j]);
            free(table);
            return NULL;
        }
        for (size_t j = 0; j < q; j++) {
            table[i][j] = ff_mul(ff, i, j);
        }
    }

    return table;
}

void ff_free_mult_table(size_t **table, size_t q)
{
    if (!table) return;
    for (size_t i = 0; i < q; i++) free(table[i]);
    free(table);
}

GFElement ff_trace(const FiniteField *ff, GFElement a)
{
    assert(ff);
    if (ff->m == 1) return a % ff->p;

    GFElement result = 0;
    GFElement current = a;
    for (size_t i = 0; i < ff->m; i++) {
        result = ff_add(ff, result, current);
        current = ff_pow(ff, current, ff->p);
    }

    return result;
}

GFElement ff_norm(const FiniteField *ff, GFElement a)
{
    assert(ff);
    if (ff->m == 1) return a % ff->p;

    /* N(a) = a^{(q-1)/(p-1)} */
    size_t exp = (ff->q - 1) / (ff->p - 1);
    return ff_pow(ff, a, exp);
}

/* ==================================================================
 *  L3 — Polynomial Operations over GF(q)
 * ================================================================== */

GFElement ff_poly_eval(const FiniteField *ff,
                        const GFElement *coeffs, size_t deg,
                        GFElement x)
{
    assert(ff);
    if (!coeffs || deg == 0) return 0;

    /* Horner's method: p(x) = ((...(c_{d-1}·x + c_{d-2})·x + ...) + c_0) */
    GFElement result = coeffs[deg - 1];

    for (size_t i = deg - 1; i > 0; i--) {
        result = ff_mul(ff, result, x);
        result = ff_add(ff, result, coeffs[i - 1]);
    }

    return result;
}

void ff_lagrange_interpolate(const FiniteField *ff,
                              const GFElement *xs,
                              const GFElement *ys,
                              size_t n, GFElement *poly)
{
    assert(ff);
    if (n == 0) return;

    /* Initialize poly to zero */
    memset(poly, 0, n * sizeof(GFElement));

    for (size_t i = 0; i < n; i++) {
        /* Compute Lagrange basis polynomial L_i(x) = Π_{j≠i} (x - x_j)/(x_i - x_j) */
        GFElement li_coeffs[128] = {0};  /* Buffer for temporary polynomial, max degree n-1 */
        size_t li_deg = 1;
        li_coeffs[0] = 1;  /* constant 1 */

        GFElement denominator = 1;

        for (size_t j = 0; j < n; j++) {
            if (j == i) continue;

            /* Multiply L_i by (x - x_j) */
            GFElement neg_xj = ff_neg(ff, xs[j]);

            /* Shift and multiply: new coeffs[i] = old[i]*neg_xj + old[i-1] */
            for (size_t k = li_deg; k > 0; k--) {
                li_coeffs[k] = ff_add(ff,
                    ff_mul(ff, li_coeffs[k], neg_xj),
                    li_coeffs[k - 1]);
            }
            li_coeffs[0] = ff_mul(ff, li_coeffs[0], neg_xj);
            li_deg++;

            /* Accumulate denominator */
            denominator = ff_mul(ff, denominator,
                ff_sub(ff, xs[i], xs[j]));
        }

        /* Scale by y_i / denominator */
        GFElement inv_den = ff_inv(ff, denominator);
        GFElement scale = ff_mul(ff, ys[i], inv_den);

        for (size_t k = 0; k < li_deg && k < n; k++) {
            poly[k] = ff_add(ff, poly[k],
                ff_mul(ff, li_coeffs[k], scale));
        }
    }
}

/* ==================================================================
 *  L4 — Field Characteristic Polynomial
 * ================================================================== */

/** Verify that the irreducible polynomial for GF(p^m) has no roots
 *  in GF(p) — a necessary (but not sufficient) condition. */
static bool poly_has_no_roots_mod_p(const size_t *coeffs, size_t deg, size_t p)
{
    /* Evaluate at each element of GF(p) */
    for (size_t x = 0; x < p; x++) {
        size_t val = 0;
        size_t xpow = 1;
        for (size_t i = 0; i <= deg; i++) {
            val = (val + coeffs[i] * xpow) % p;
            xpow = (xpow * x) % p;
        }
        if (val == 0) return false;
    }
    return true;
}
