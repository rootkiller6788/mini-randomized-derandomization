/*
 * list_decode_core.h — Core Definitions for List-Decodable Codes
 *
 * List-decodable codes (LDC) generalise unique decoding by allowing
 * the decoder to output a small list of candidate codewords. This
 * paradigm, pioneered by Elias (1957) and Wozencraft (1958),
 * was revived by Sudan (1997) and Guruswami-Sudan (1999).
 *
 * References:
 *   - Guruswami, "List Decoding of Error-Correcting Codes" (2004)
 *   - Sudan, J. Complexity 13(1):180-193, 1997
 *   - Johnson, IEEE Trans. Info. Theory 8(5):29-37, 1962
 */
#ifndef LD_CORE_H
#define LD_CORE_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

/* ------------------------------------------------------------------ */
/*  L1 — Core Type Definitions                                        */
/* ------------------------------------------------------------------ */

/** Generic code parameters: [n,k,d]_q where
 *  n = block length, k = message length, d = minimum distance,
 *  q = alphabet size (must be a prime power for linear codes). */
typedef struct {
    size_t n;      /* block length */
    size_t k;      /* message/information length */
    size_t d;      /* minimum distance */
    size_t q;      /* alphabet size (prime power for GF(q)) */
} CodeParams;

/** A single codeword — a sequence of alphabet symbols. */
typedef struct {
    size_t *symbols;   /* array of length `len` */
    size_t  len;       /* number of symbols */
    size_t  alphabet;  /* q — alphabet size */
} Codeword;

/** A (possibly empty) list of codewords returned by a list decoder. */
typedef struct {
    Codeword *words;   /* dynamically allocated array */
    size_t    count;   /* number of codewords in list */
} CodewordList;

/** Parameters controlling a single list-decoding invocation.
 *  `delta` = fraction of errors tolerated (0 ≤ δ ≤ 1).
 *  `list_size` = an upper bound on the output list size. */
typedef struct {
    size_t n;           /* block length */
    double delta;       /* error fraction */
    size_t list_size;   /* max number of candidate codewords */
} ListDecodeParams;

/** Euclidean / Hamming distance result between two codewords. */
typedef struct {
    size_t hamming;     /* plain Hamming distance */
    double relative;    /* d_H / n */
} CodeDistance;

/* ------------------------------------------------------------------ */
/*  L2 — Core Concepts: Encoding, Decoding, Bounds                    */
/* ------------------------------------------------------------------ */

/**
 * Encode a binary message into a codeword over the code's alphabet.
 * @param cp   Code parameters [n,k,d,q]
 * @param msg  Binary message bits (length k·ceil(log2 q) or k symbols)
 * @param msg_len Number of message symbols
 * @return A newly allocated Codeword of length n over alphabet q.
 * Complexity: O(n·k) time.
 */
Codeword ld_encode(const CodeParams *cp, const size_t *msg, size_t msg_len);

/**
 * List-decode a received word: return *all* codewords within
 * relative Hamming distance δ = lp->delta of `received`.
 * @param received  Received word of length n
 * @param lp        List-decode parameters (n, δ, list_size)
 * @param max_list  Hard upper bound on output list cardinality
 * @return A CodewordList containing 0..max_list candidates.
 * Complexity: O(n·|C|) brute force (theoretical reference).
 */
CodewordList ld_list_decode(const Codeword *received,
                             const ListDecodeParams *lp,
                             size_t max_list);

/* ------------------------------------------------------------------ */
/*  L4 — Fundamental Bounds (C implementation)                        */
/* ------------------------------------------------------------------ */

/**
 * Johnson bound for list decoding.
 * For a code of block length n and minimum distance d over alphabet
 * of size q, the Johnson radius e = n·(1 - sqrt(1 - d/n)) gives the
 * maximum error radius for which the list size is polynomially bounded.
 *
 * Theorem (Johnson 1962): Any code C ⊆ [q]^n with distance d has
 *   |B(x, e) ∩ C| ≤ q·n   for e ≤ n·(1 - sqrt{1 - d/n}).
 *
 * @return The Johnson bound list size L_J(n,d,q).
 */
size_t ld_johnson_bound(size_t n, size_t d, size_t q);

/**
 * List-decoding radius à la Johnson.
 * @return floor(n·(1 - sqrt(1 - (double)d/n))).
 */
double ld_johnson_radius(size_t n, size_t d);

/**
 * Combinatorial list-decoding capacity (asymptotic).
 * For an (n,k) code, a fraction δ of errors can be list-decoded
 * with polynomial list size iff δ ≤ 1 - R - ε.
 * @return The capacity limit error fraction for rate R = k/n.
 */
double ld_list_decoding_capacity(double rate, double epsilon);

/**
 * Griesmer bound: for a linear [n,k,d]_q code,
 *   n ≥ Σ_{i=0}^{k-1} ⌈ d / q^i ⌉.
 * @return Lower bound on n given k, d, q.
 */
size_t ld_griesmer_bound(size_t k, size_t d, size_t q);

/* ------------------------------------------------------------------ */
/*  L2 (continued) — Distance & Validity                             */
/* ------------------------------------------------------------------ */

/** Compute Hamming distance (absolute and relative) between two codewords. */
CodeDistance ld_hamming_distance(const Codeword *a, const Codeword *b);

/** Agnostically test whether a code is list-decodable up to radius r
 *  with list size at most L (exhaustive check for small codes). */
bool ld_is_list_decodable(const CodeParams *cp, double radius, size_t L);

/** Free a singly-allocated codeword. */
void ld_free_codeword(Codeword *cw);

/** Free an entire list of codewords. */
void ld_free_list(CodewordList *cl);

/** Deep copy a codeword. */
Codeword ld_copy_codeword(const Codeword *src);

/** Print a codeword to stdout (for debugging). */
void ld_print_codeword(const Codeword *cw);

#endif /* LD_CORE_H */
