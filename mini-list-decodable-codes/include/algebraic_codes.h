/*
 * algebraic_codes.h — Algebraic Code Structures
 *
 * Defines linear codes, cyclic codes, BCH codes, and LDPC codes
 * — the algebraic backbone for many list-decoding constructions.
 *
 * Linear codes over GF(q): a k-dimensional subspace of GF(q)^n.
 * Represented by a generator matrix G (k×n) or parity-check matrix H.
 *
 * References:
 *   - MacWilliams & Sloane, "The Theory of Error-Correcting Codes" (1977)
 *   - Gallager, "Low-Density Parity-Check Codes" (1963)
 *   - Hocquenghem (1959), Bose & Ray-Chaudhuri (1960) for BCH
 */
#ifndef ALG_CODES_H
#define ALG_CODES_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

/* ------------------------------------------------------------------ */
/*  L1 — Linear Code [n,k,d]_q                                       */
/* ------------------------------------------------------------------ */

/** A linear code is a k-dimensional subspace of GF(q)^n.
 *  generator: k×n matrix over GF(q). */
typedef struct {
    size_t  n;             /* block length */
    size_t  k;             /* dimension */
    size_t  q;             /* field size (prime power) */
    size_t  d;             /* computed minimum distance */
    size_t **generator;    /* k×n generator matrix (row-major) */
} LinearCode;

/** A cyclic code of length n over GF(q).
 *  Defined by a generator polynomial g(x) that divides x^n - 1.
 *  dimension k = n - deg(g). */
typedef struct {
    size_t  n;                /* block length */
    size_t  k;                /* dimension = n - deg(generator_poly) */
    size_t  q;                /* field size */
    size_t  deg_gen;          /* degree of generator polynomial */
    size_t *generator_poly;   /* coefficients of g(x), degree d */
} CyclicCode;

/** Low-Density Parity-Check (LDPC) code.
 *  Defined by a sparse (r × n) parity-check matrix H.
 *  Regular LDPC: each column has weight wc, each row weight wr. */
typedef struct {
    size_t  n;              /* block length */
    size_t  k;              /* approximate dimension ≈ n - r (if H full rank) */
    size_t  r;              /* number of parity-check equations (rows of H) */
    size_t  wc;             /* column weight (regular) */
    size_t  wr;             /* row weight (regular) */
    size_t **parity_check;  /* r × n sparse matrix over GF(2) */
} LDPCCode;

/* ------------------------------------------------------------------ */
/*  L2 — Linear Code Operations                                       */
/* ------------------------------------------------------------------ */

/** Create a blank linear code (allocates generator/parity arrays).
 *  Caller fills generator matrix row by row. */
LinearCode *lc_create(size_t n, size_t k, size_t q);

/** Create a linear code from a given generator matrix G (k×n, row-major). */
LinearCode *lc_from_generator(size_t n, size_t k, size_t q,
                               const size_t *G_flat);

/** Encode message msg (k symbols) by multiplying by generator: c = msg × G. */
size_t lc_encode(const LinearCode *lc, const size_t *msg, size_t msg_len,
                 size_t *codeword_out);

/** Bounded-distance decoder: find codeword within Hamming distance t. */
bool lc_decode_bounded_distance(const LinearCode *lc,
                                const size_t *received,
                                size_t t, size_t *decoded_out);

/** List-decode a linear code: return all codewords within radius pn. */
size_t lc_list_decode(const LinearCode *lc, const size_t *received,
                      double radius, size_t **decoded_out,
                      size_t max_list);

/** Compute minimum distance of a linear code (brute force over k small). */
size_t lc_compute_distance(const LinearCode *lc);

/** Free a linear code. */
void lc_free(LinearCode *lc);

/* ------------------------------------------------------------------ */
/*  L2 — Cyclic Code Operations                                       */
/* ------------------------------------------------------------------ */

/** Create a BCH code with designed distance δ=2t+1 correcting t errors.
 *  Uses a primitive root α to build the generator polynomial. */
CyclicCode *cc_create_bch(size_t n, size_t t, size_t q);

/** Encode with a cyclic code (polynomial multiplication). */
size_t cc_encode(const CyclicCode *cc, const size_t *msg, size_t msg_len,
                 size_t *codeword_out);

/** List-decode a cyclic code using Sudan-type polynomial reconstruction. */
size_t cc_list_decode(const CyclicCode *cc, const size_t *received,
                      double radius, size_t **decoded_out,
                      size_t max_list);

/** Free a cyclic code. */
void cc_free(CyclicCode *cc);

/* ------------------------------------------------------------------ */
/*  L5 — LDPC Code Operations                                         */
/* ------------------------------------------------------------------ */

/** Create an (n, wc, wr)-regular LDPC code with a random-like Tanner graph. */
LDPCCode *ldpc_create_regular(size_t n, size_t wc, size_t wr);

/** Belief-propagation (sum-product) decoder for LDPC codes.
 *  @param llr       Log-likelihood ratios (length n).
 *  @param decoded   Output hard decisions (length n).
 *  @param max_iter  Maximum number of message-passing iterations.
 *  @return true if decoding converged to a valid codeword. */
bool ldpc_belief_propagation_decode(const LDPCCode *ldpc,
                                    const double *llr,
                                    size_t *decoded,
                                    size_t max_iter);

/** Free an LDPC code. */
void ldpc_free(LDPCCode *ldpc);

#endif /* ALG_CODES_H */
