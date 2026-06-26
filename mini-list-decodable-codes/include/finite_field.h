/*
 * finite_field.h — Finite Field GF(q) Arithmetic
 *
 * Implements arithmetic in prime fields GF(p) and prime-power
 * fields GF(p^m).  Essential for all algebraic coding theory.
 *
 * For simplicity and portability (C99), the primary implementation
 * is for prime fields GF(p) using modular arithmetic.
 * Extension fields GF(p^m) use polynomial representation.
 *
 * References:
 *   - Lidl & Niederreiter, "Finite Fields" (2nd ed., 1997)
 *   - McEliece, "Finite Fields for Computer Scientists and Engineers" (1987)
 */
#ifndef FINITE_FIELD_H
#define FINITE_FIELD_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>

/* ------------------------------------------------------------------ */
/*  L1 / L3 — Field Element and Field Structure                       */
/* ------------------------------------------------------------------ */

/** A finite field element.  For prime fields it's an integer mod p.
 *  For extension fields GF(p^m) it can represent polynomial coefficients. */
typedef size_t GFElement;

/** Finite field descriptor.
 *  p      = prime characteristic
 *  m      = extension degree (1 = prime field)
 *  q      = p^m (total number of elements)
 *  irreducible_poly = for m > 1, coefficients of the irreducible polynomial
 *                     over GF(p) of degree m, used for multiplication.
 */
typedef struct {
    size_t  p;                    /* prime characteristic */
    size_t  m;                    /* extension degree */
    size_t  q;                    /* p^m = field size */
    size_t *irreducible_poly;     /* coefficients for GF(p^m) when m > 1 */
} FiniteField;

/* ------------------------------------------------------------------ */
/*  L1 — Field Creation and Destruction                               */
/* ------------------------------------------------------------------ */

/** Create a prime field GF(p).
 *  p must be prime.  (Primality is checked.) */
FiniteField *ff_create_prime(size_t p);

/** Create an extension field GF(p^m) given an irreducible polynomial.
 *  irr_poly has m+1 coefficients: irr_poly[0] + irr_poly[1]·x + ... */
FiniteField *ff_create_extension(size_t p, size_t m,
                                  const size_t *irr_poly);

/** Free a finite field. */
void ff_free(FiniteField *ff);

/* ------------------------------------------------------------------ */
/*  L2 — Field Arithmetic                                            */
/* ------------------------------------------------------------------ */

/** Add a + b in GF(q). */
GFElement ff_add(const FiniteField *ff, GFElement a, GFElement b);

/** Subtract a - b in GF(q). */
GFElement ff_sub(const FiniteField *ff, GFElement a, GFElement b);

/** Multiply a × b in GF(q). */
GFElement ff_mul(const FiniteField *ff, GFElement a, GFElement b);

/** Compute a^{-1} in GF(q).  Returns 0 if a == 0 (no inverse). */
GFElement ff_inv(const FiniteField *ff, GFElement a);

/** Divide a / b in GF(q).  Returns 0 if b == 0. */
GFElement ff_div(const FiniteField *ff, GFElement a, GFElement b);

/** Compute a^e mod q efficiently (square-and-multiply). */
GFElement ff_pow(const FiniteField *ff, GFElement a, size_t exp);

/** Test equality a == b in GF(q). */
bool ff_eq(const FiniteField *ff, GFElement a, GFElement b);

/** Negate -a in GF(q). */
GFElement ff_neg(const FiniteField *ff, GFElement a);

/* ------------------------------------------------------------------ */
/*  L2 — Utility Functions                                            */
/* ------------------------------------------------------------------ */

/** Check if p is prime (simple trial division, for small p). */
bool ff_is_prime(size_t p);

/** Find a primitive element (generator) of the multiplicative group GF(q)^*.
 *  Uses a simple search with factorization of q-1. */
GFElement ff_primitive_element(const FiniteField *ff);

/** Generate the multiplicative table for small fields (a × b lookup).
 *  Useful for speeding up repeated operations. */
size_t **ff_build_mult_table(const FiniteField *ff);

/** Free a multiplication table. */
void ff_free_mult_table(size_t **table, size_t q);

/** Compute the trace of an element: Tr(a) = a + a^p + a^{p^2} + ... + a^{p^{m-1}}. */
GFElement ff_trace(const FiniteField *ff, GFElement a);

/** Compute the norm: N(a) = a^{(q-1)/(p-1)}. */
GFElement ff_norm(const FiniteField *ff, GFElement a);

/* ------------------------------------------------------------------ */
/*  L3 — Polynomial Operations over GF(q)                             */
/* ------------------------------------------------------------------ */

/**
 * Evaluate polynomial p(x) = p_0 + p_1·x + ... + p_{deg-1}·x^{deg-1}
 * at point x ∈ GF(q).
 * @param coeffs  Coefficient array, length deg.
 * @param deg     Polynomial degree (number of coefficients).
 * @param x       Evaluation point.
 * @return p(x).
 */
GFElement ff_poly_eval(const FiniteField *ff,
                        const GFElement *coeffs, size_t deg,
                        GFElement x);

/**
 * Lagrange interpolation: given points (x_i, y_i) for i=0..n-1,
 * find the unique polynomial of degree < n that passes through all points.
 * @param xs   n evaluation points (distinct).
 * @param ys   n values.
 * @param n    Number of points.
 * @param poly Output polynomial coefficients (length n), caller allocates.
 */
void ff_lagrange_interpolate(const FiniteField *ff,
                              const GFElement *xs,
                              const GFElement *ys,
                              size_t n, GFElement *poly);

#endif /* FINITE_FIELD_H */