/*
 * polynomial.h — Univariate & Bivariate Polynomial Operations
 *
 * Polynomial manipulation over finite fields is the algorithmic
 * core of algebraic (list-)decoding.  This module provides:
 *   - Univariate polynomial arithmetic (+, -, ×, ÷, GCD)
 *   - Bivariate polynomial representation for Sudan-type algorithms
 *   - Root-finding over finite fields (Chien search)
 *   - Hasse derivative for multiplicity-based decoding
 *
 * References:
 *   - von zur Gathen & Gerhard, "Modern Computer Algebra" (3rd ed., 2013)
 *   - Sudan, J. Complexity 13(1):180-193, 1997
 *   - Guruswami & Sudan, IEEE Trans. Info. Theory 45(6):1757-1767, 1999
 */
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include "finite_field.h"

/* ------------------------------------------------------------------ */
/*  L3 — Univariate Polynomial Types                                  */
/* ------------------------------------------------------------------ */

/** A univariate polynomial over GF(q): p(x) = c_0 + c_1·x + ... + c_{d-1}·x^{d-1}.
 *  degree = highest index with non-zero coefficient, or 0 for zero poly.
 *  coeffs[i] = coefficient of x^i. */
typedef struct {
    GFElement *coeffs;   /* coefficients[0..deg] */
    size_t     deg;      /* degree (coeffs[deg] ≠ 0 for non-zero poly) */
    size_t     capacity; /* allocated size (≥ deg+1) */
} UnivariatePoly;

/* ------------------------------------------------------------------ */
/*  L2 — Univariate Polynomial Operations                              */
/* ------------------------------------------------------------------ */

/** Create zero polynomial "0". */
UnivariatePoly *upoly_zero(void);

/** Create a constant polynomial "c". */
UnivariatePoly *upoly_constant(const FiniteField *ff, GFElement c);

/** Create polynomial from coefficient array.  Takes ownership of coeffs. */
UnivariatePoly *upoly_from_coeffs(GFElement *coeffs, size_t deg, size_t capacity);

/** Deep copy a polynomial. */
UnivariatePoly *upoly_copy(const UnivariatePoly *p);

/** Free a polynomial. */
void upoly_free(UnivariatePoly *p);

/** Trim trailing zero coefficients and update degree. */
void upoly_trim(UnivariatePoly *p);

/** Polynomial addition: r = a + b. */
UnivariatePoly *upoly_add(const FiniteField *ff,
                           const UnivariatePoly *a,
                           const UnivariatePoly *b);

/** Polynomial subtraction: r = a - b. */
UnivariatePoly *upoly_sub(const FiniteField *ff,
                           const UnivariatePoly *a,
                           const UnivariatePoly *b);

/** Polynomial multiplication: r = a × b.  O(deg(a)·deg(b)). */
UnivariatePoly *upoly_mul(const FiniteField *ff,
                           const UnivariatePoly *a,
                           const UnivariatePoly *b);

/** Polynomial evaluation: p(x) at point x. */
GFElement upoly_eval(const FiniteField *ff,
                      const UnivariatePoly *p, GFElement x);

/** Formal derivative p'(x). */
UnivariatePoly *upoly_derivative(const FiniteField *ff,
                                  const UnivariatePoly *p);

/** Extended Euclidean algorithm: returns gcd(a,b) and sets *s, *t
 *  such that s·a + t·b = gcd(a,b). */
UnivariatePoly *upoly_extended_gcd(const FiniteField *ff,
                                    const UnivariatePoly *a,
                                    const UnivariatePoly *b,
                                    UnivariatePoly **s,
                                    UnivariatePoly **t);

/** Polynomial division: q = a / b, r = a mod b. */
void upoly_divmod(const FiniteField *ff,
                  const UnivariatePoly *a,
                  const UnivariatePoly *b,
                  UnivariatePoly **q,
                  UnivariatePoly **r);

/** Find all roots of p(x) in GF(q) by Chien search.
 *  @param roots_out  Output: array of roots (caller frees).
 *  @return Number of roots found. */
size_t upoly_find_roots(const FiniteField *ff,
                         const UnivariatePoly *p,
                         GFElement **roots_out);

/* ------------------------------------------------------------------ */
/*  L5 — Bivariate Polynomials for Sudan-type Algorithms               */
/* ------------------------------------------------------------------ */

/**
 * A bivariate polynomial Q(x,y) = Σ a_{i,j}·x^i·y^j.
 * Stored as a flat array of size (deg_x+1)×(deg_y+1) in row-major order.
 * coeffs[i*(deg_y+1) + j] = coefficient of x^i·y^j.
 */
typedef struct {
    GFElement *coeffs;   /* flat 2D array */
    size_t     deg_x;    /* max degree in x */
    size_t     deg_y;    /* max degree in y */
} BivariatePoly;

/** Create a zero bivariate polynomial. */
BivariatePoly *bpoly_create(size_t deg_x, size_t deg_y);

/** Free a bivariate polynomial. */
void bpoly_free(BivariatePoly *bp);

/** Get coefficient of x^i·y^j. */
GFElement bpoly_coeff(const BivariatePoly *bp, size_t i, size_t j);

/** Set coefficient of x^i·y^j. */
void bpoly_set_coeff(BivariatePoly *bp, size_t i, size_t j, GFElement val);

/**
 * Evaluate bivariate polynomial at (x0, y0).
 * Q(x0, y0) = Σ_{i,j} a_{i,j}·x0^i·y0^j.
 */
GFElement bpoly_eval(const FiniteField *ff,
                      const BivariatePoly *bp,
                      GFElement x0, GFElement y0);

/**
 * Hasse derivative of Q(x,y) with respect to x (order a) and y (order b).
 * D_x^a D_y^b Q evaluated at (x0, y0).
 */
GFElement bpoly_hasse_derivative(const FiniteField *ff,
                                  const BivariatePoly *bp,
                                  size_t a, size_t b,
                                  GFElement x0, GFElement y0);

/**
 * Factor bivariate polynomial Q(x,y) over the univariate factors y - p(x).
 * This is the core of Sudan's algorithm: find all polynomials p
 * such that (y - p(x)) divides Q(x,y).
 *
 * @param factors_out  Array of UnivariatePoly* (caller frees each).
 * @return Number of factors found.
 */
size_t bpoly_factor_linear_y(const FiniteField *ff,
                              const BivariatePoly *Q,
                              UnivariatePoly ***factors_out);

/* ------------------------------------------------------------------ */
/*  L3 — Monomial Ordering & Weighted Degree                           */
/* ------------------------------------------------------------------ */

/** (1,k-1)-weighted degree of monomial x^i·y^j:
 *  wdeg(i,j) = i + j·(k-1) for RS codes with message length k. */
size_t bpoly_weighted_deg(size_t i, size_t j, size_t k);

/** Total (1,k-1)-weighted degree of bivariate polynomial. */
size_t bpoly_total_weighted_deg(const BivariatePoly *bp, size_t k);

#endif /* POLYNOMIAL_H */