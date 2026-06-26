/*
 * rs_codes.h — Reed-Solomon Codes and Their List-Decoding Algorithms
 *
 * Reed-Solomon (RS) codes are the most important family of algebraic
 * error-correcting codes.  For any n ≤ q and k ≤ n, the [n,k,n-k+1]_q
 * RS code evaluates polynomials of degree < k at n distinct points.
 *
 * Key results for list decoding:
 *   - Sudan (1997): poly-time list decoding beyond d/2
 *   - Guruswami-Sudan (1999): multiplicity-based, better radius
 *   - Guruswami-Rudra (2008): folded RS codes achieve capacity
 *
 * References:
 *   - Reed & Solomon, J. SIAM 8(2):300-304, 1960
 *   - Sudan, J. Complexity 13(1):180-193, 1997
 *   - Guruswami & Sudan, IEEE Trans. Info. Theory 45(6):1757-1767, 1999
 */
#ifndef RS_CODES_H
#define RS_CODES_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

/* ------------------------------------------------------------------ */
/*  L1 — Reed-Solomon Code Structure                                  */
/* ------------------------------------------------------------------ */

/** Classical Reed-Solomon code over GF(q) with evaluation points.
 *  n = block length (≤ q), k = message length, q = field size (prime power).
 *  evaluation_points[i] ∈ GF(q) are distinct. */
typedef struct {
    size_t  n;                  /* block length */
    size_t  k;                  /* message (information) symbols */
    size_t  q;                  /* field size (must be a prime power) */
    size_t *evaluation_points;  /* array[0..n-1] of distinct field elements */
} RSCode;

/* ------------------------------------------------------------------ */
/*  L1 — Reed-Solomon Code Operations                                 */
/* ------------------------------------------------------------------ */

/** Create an [n,k,n-k+1]_q RS code with distinct evaluation points.
 *  If eval_pts == NULL, uses standard points {0,1,2,...,n-1} mod q. */
RSCode *rs_create(size_t n, size_t k, size_t q, const size_t *eval_pts);

/** Encode a message (k coefficients of a polynomial of degree < k)
 *  into an n-symbol codeword by evaluating the polynomial at all
 *  evaluation points.  codeword[i] = f(eval_pt[i]) mod q.
 *  @return Number of symbols written (n). */
size_t rs_encode(const RSCode *rs, const size_t *msg, size_t msg_len,
                 size_t *codeword_out);

/** Compute the minimum distance of the code = n - k + 1. */
size_t rs_min_distance(const RSCode *rs);

/** Free an RS code structure. */
void rs_free(RSCode *rs);

/* ------------------------------------------------------------------ */
/*  L5 — List-Decoding Algorithms for RS Codes                        */
/* ------------------------------------------------------------------ */

/**
 * Sudan's list-decoding algorithm (1997).
 * Given a received word (r_0,...,r_{n-1}), find *all* polynomials
 * f of degree < k such that f(eval_pt[i]) = r_i for at least
 * T = ceil(sqrt(2·(k-1)·n)) distinct positions i.
 *
 * @param multiplicity  Controls the interpolation step (≥ 1).
 * @param decoded_out   Output: array of size L*k, each block of k symbols
 *                      is one candidate polynomial's coefficients.
 * @return Number L of candidate polynomials found.
 */
size_t rs_sudan_list_decode(const RSCode *rs, const size_t *received,
                            size_t multiplicity,
                            size_t **decoded_out);

/**
 * Guruswami-Sudan (1999) multiplicity-based list decoder.
 * Extends Sudan by allowing each interpolation constraint to
 * have multiplicity s ≥ 1, achieving a larger decoding radius.
 *
 * @param s  multiplicity parameter (≥ 1, higher = better radius, slower).
 * @param decoded_out   Output: array of size L*k.
 * @return Number L of candidate polynomials.
 */
size_t rs_guruswami_sudan(const RSCode *rs, const size_t *received,
                          size_t s,
                          size_t **decoded_out);

/**
 * Welch-Berlekamp unique decoder (reference baseline).
 * Corrects up to floor((n-k)/2) errors using rational interpolation.
 * @param decoded_out  Array of k symbols (single polynomial).
 * @return true on success, false if too many errors.
 */
bool rs_welch_berlekamp(const RSCode *rs, const size_t *received,
                        size_t *decoded_out);

/* ------------------------------------------------------------------ */
/*  L4 — List-Decoding Bounds for RS Codes                            */
/* ------------------------------------------------------------------ */

/**
 * Johnson radius for RS codes:
 *   e_J = n·(1 - sqrt(k/n)) = n·(1 - sqrt{R}).
 * The maximum error count for which list size is polynomially bounded.
 */
double rs_johnson_radius(const RSCode *rs);

/**
 * Johnson bound for RS codes:
 * List size ≤ q·n for errors ≤ Johnson radius.
 * @return Upper bound on list size for given number of errors.
 */
size_t rs_johnson_bound_list_size(const RSCode *rs, size_t errors);

/**
 * List-decoding radius for RS codes under Sudan's algorithm:
 *   e_S(n,k) = n - sqrt(2·(k-1)·n).
 */
double rs_sudan_radius(const RSCode *rs);

/**
 * Guruswami-Sudan decoding radius for multiplicity s:
 *   e_GS(n,k,s) ≈ n·(1 - sqrt{(k-1)/n}·sqrt{1+1/s}).
 */
double rs_gs_radius(const RSCode *rs, size_t multiplicity);

/**
 * Capacity limit for list-decoding RS codes over GF(q):
 *   δ_cap = 1 - R.
 * For δ < δ_cap, exp(O(1/ε)) list size is asymptotically possible.
 */
double rs_list_capacity_limit(size_t q, double rate);

/* ------------------------------------------------------------------ */
/*  L2 — Helper: Verify List-Decoding Properties                      */
/* ------------------------------------------------------------------ */

/**
 * Verify the Johnson bound for a specific RS instance.
 * Checks: for all codewords c, |B(r, e) ∩ C| ≤ L_J.
 * @return true if bound holds.
 */
bool rs_verify_johnson_bound(const RSCode *rs, const size_t *received,
                             size_t errors, size_t list_cap);

#endif /* RS_CODES_H */
