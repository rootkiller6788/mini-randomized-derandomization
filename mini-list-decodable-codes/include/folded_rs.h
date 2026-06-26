/*
 * folded_rs.h — Folded Reed-Solomon Codes
 *
 * Folded Reed-Solomon (FRS) codes, introduced by Guruswami & Rudra
 * (2008), are the first explicit construction of codes that achieve
 * the list-decoding capacity.
 *
 * Construction: take an RS code over GF(q) with block length n,
 * then "fold" consecutive symbols into bundles of size m, yielding
 * a new code over alphabet GF(q)^m with block length N = n/m.
 *
 * Key Theorem (Guruswami-Rudra 2008):
 *   For any ε > 0, there exist explicit folded RS codes of rate R
 *   that are list-decodable up to error fraction 1 - R - ε
 *   with list size poly(n, 1/ε).
 *
 * References:
 *   - Guruswami & Rudra, IEEE Trans. Info. Theory 54(2):692-706, 2008
 *   - Guruswami & Xing, "List Decoding of Folded RS Codes via Linear-Algebraic
 *     Approach" (STOC 2012)
 *   - Dvir & Lovett, "Subspace Evasive Sets" (STOC 2012)
 */
#ifndef FOLDED_RS_H
#define FOLDED_RS_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

/* ------------------------------------------------------------------ */
/*  L1 — Folded Reed-Solomon Code Structure                            */
/* ------------------------------------------------------------------ */

/**
 * A folded RS code: [N, K]_Q where
 *   N = n / m  (folded block length),
 *   K = k / m  (folded message length),
 *   Q = q^m    (folded alphabet size).
 *
 * The original (unfolded) RS code has parameters [n, k, n-k+1]_q.
 */
typedef struct {
    size_t n;                /* unfolded block length */
    size_t k;                /* unfolded message length */
    size_t q;                /* base field size */
    size_t m;                /* folding parameter (bundle size) */
    size_t N;                /* folded block length = n / m */
    size_t K;                /* folded message length = k / m */
    size_t *eval_points;     /* evaluation points (length n) */
} FoldedRSCode;

/* ------------------------------------------------------------------ */
/*  L2 — FRS Code Operations                                          */
/* ------------------------------------------------------------------ */

/**
 * Create a folded RS code.
 * @param n  Unfolded block length (must be divisible by m).
 * @param k  Unfolded message length (must be divisible by m).
 * @param q  Base field size GF(q).
 * @param m  Folding parameter.
 * @param eval_pts  Evaluation points (length n, or NULL for standard).
 */
FoldedRSCode *frs_create(size_t n, size_t k, size_t q, size_t m,
                          const size_t *eval_pts);

/** Encode a message into a folded codeword.
 *  msg has k symbols over GF(q).  Output codeword has N blocks
 *  of m symbols each (N·m = n total symbols). */
void frs_encode(const FoldedRSCode *frs, const size_t *msg,
                size_t **codeword_out);

/** Free a folded RS code. */
void frs_free(FoldedRSCode *frs);

/* ------------------------------------------------------------------ */
/*  L5 — List-Decoding Algorithm for Folded RS Codes                   */
/* ------------------------------------------------------------------ */

/**
 * List-decode a folded RS code to capacity.
 *
 * Uses the linear-algebraic approach (Guruswami-Xing 2012):
 * reduce to solving a linear system over the base field.
 *
 * @param received     Received word (N·m symbols over GF(q)).
 * @param epsilon      Target error margin δ = 1 - R - ε.
 * @param decoded_out  Output: array of candidate messages (each of length k).
 * @param max_list     Maximum number of candidates.
 * @return Number of candidates found.
 */
size_t frs_list_decode_to_capacity(const FoldedRSCode *frs,
                                    const size_t *received,
                                    double epsilon,
                                    size_t **decoded_out,
                                    size_t max_list);

/**
 * Compute the list-decoding radius for folded RS.
 *   δ_max(m) = (m/(m+1))·(1 - k/n·(m+1)/(m-R)).
 * For m → ∞, δ_max → 1 - R (capacity).
 */
double frs_list_decoding_radius(const FoldedRSCode *frs);

/**
 * Verify that the code achieves the claimed error correction.
 * @return true if the list-decoding radius is within ε of capacity.
 */
bool frs_verify_capacity_approaching(const FoldedRSCode *frs, double epsilon);

/* ------------------------------------------------------------------ */
/*  L4 — Theoretical Bounds for Folded RS                              */
/* ------------------------------------------------------------------ */

/**
 * List-size bound for folded RS codes.
 * Theorem (Guruswami-Rudra): For δ < 1-R-ε, list size ≤ (n/ε)^{O(1/ε)}.
 */
size_t frs_list_size_bound(const FoldedRSCode *frs, double epsilon);

/** Minimum distance of the folded RS code. */
size_t frs_min_distance(const FoldedRSCode *frs);

/** Rate of the folded code: R_folded = K/N = k/n. */
double frs_rate(const FoldedRSCode *frs);

/** Test whether given parameters can be list-decoded up to capacity. */
bool frs_is_capacity_achieving(size_t n, size_t k, size_t m, double target_delta);

#endif /* FOLDED_RS_H */