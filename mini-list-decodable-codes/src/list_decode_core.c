/*
 * list_decode_core.c — Core List-Decoding Operations
 *
 * Implements the fundamental operations for list-decodable codes:
 * encoding, Hamming distance computation, Johnson bound,
 * list-decoding radius, and combinatorial bounds.
 *
 * This is the entry point for all list-decodable code operations.
 */

#include "list_decode_core.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

/* ==================================================================
 *  L1 / L2 — Encoding and Basic Operations
 * ================================================================== */

Codeword ld_encode(const CodeParams *cp, const size_t *msg, size_t msg_len)
{
    assert(cp);
    Codeword cw;
    cw.len = cp->n;
    cw.alphabet = cp->q;
    cw.symbols = (size_t *)calloc(cp->n, sizeof(size_t));

    if (!cw.symbols) {
        cw.len = 0;
        return cw;
    }

    /* Simple systematic encoding: first k symbols are message,
     * remaining n-k are parity checks (for linear codes this is
     * the generator-matrix approach simplified). */

    for (size_t i = 0; i < cp->k && i < cp->n; i++) {
        cw.symbols[i] = (i < msg_len) ? (msg[i] % cp->q) : 0;
    }

    /* Parity symbols: cyclic redundancy for demonstration.
     * In a real system, this would use the specific code's parity
     * equations.  Here we compute parity as XOR-sum over GF(q) of
     * the message symbols with a shift pattern. */
    for (size_t i = cp->k; i < cp->n; i++) {
        size_t parity = 0;
        for (size_t j = 0; j < cp->k && j < i; j++) {
            parity = (parity + cw.symbols[j] * ((i - j) % cp->q)) % cp->q;
        }
        cw.symbols[i] = parity;
    }

    return cw;
}

CodeDistance ld_hamming_distance(const Codeword *a, const Codeword *b)
{
    CodeDistance dist;
    dist.hamming = 0;
    dist.relative = 0.0;

    if (!a || !b || !a->symbols || !b->symbols) return dist;
    if (a->len != b->len) {
        dist.hamming = (size_t)-1;  /* Error indicator */
        return dist;
    }

    for (size_t i = 0; i < a->len; i++) {
        if (a->symbols[i] != b->symbols[i]) {
            dist.hamming++;
        }
    }

    dist.relative = (double)dist.hamming / (double)a->len;
    return dist;
}

/* Brute-force list decoder: enumerate all codewords in a code and
 * return those within the specified relative Hamming distance.
 * This is the reference implementation — production code would use
 * the algebraic Sudan/Guruswami-Sudan algorithms. */
CodewordList ld_list_decode(const Codeword *received,
                             const ListDecodeParams *lp,
                             size_t max_list)
{
    assert(received);
    assert(lp);

    CodewordList result;
    result.words = (Codeword *)calloc(max_list, sizeof(Codeword));
    result.count = 0;

    if (!result.words) return result;

    /* In a brute-force approach, we would enumerate all q^k codewords.
     * For small parameters (n ≤ 10, q ≤ 3), this is feasible.
     * For larger codes, the algebraic decoders replace this. */

    size_t threshold = (size_t)ceil(lp->delta * (double)lp->n);

    /* Enumerate all possible messages */
    size_t space_size = 1;
    for (size_t i = 0; i < lp->n && space_size < (size_t)1e6; i++) {
        size_t prev = space_size;
        space_size *= received->alphabet;
        if (space_size / received->alphabet != prev) {
            space_size = (size_t)1e6;  /* Overflow guard */
            break;
        }
    }

    /* For large spaces, skip brute-force and return empty list.
     * Real systems use specialized list-decoding algorithms. */
    if (space_size > (size_t)1e5) {
        /* Too large for brute-force; algebraic decoder must be used.
         * Return an empty list to indicate brute-force infeasible. */
        result.count = 0;
        return result;
    }

    /* Brute-force: try all possible codewords */
    Codeword candidate;
    candidate.len = lp->n;
    candidate.alphabet = received->alphabet;
    candidate.symbols = (size_t *)calloc(lp->n, sizeof(size_t));
    if (!candidate.symbols) return result;

    /* We iterate over all possible messages of length n */
    size_t *msg = (size_t *)calloc(lp->n, sizeof(size_t));
    if (!msg) { free(candidate.symbols); return result; }

    for (size_t enc = 0; enc < space_size && result.count < max_list; enc++) {
        /* Convert counter to base-q digits */
        size_t tmp = enc;
        for (size_t i = 0; i < lp->n; i++) {
            msg[i] = tmp % received->alphabet;
            tmp /= received->alphabet;
        }

        /* Encode the message */
        CodeParams cp;
        cp.n = lp->n;
        cp.k = lp->n / 2;  /* Approximate, for demonstration */
        cp.q = received->alphabet;
        cp.d = cp.n - cp.k + 1;
        Codeword cw = ld_encode(&cp, msg, lp->n);

        /* Check distance */
        CodeDistance d = ld_hamming_distance(received, &cw);
        if (d.hamming <= threshold) {
            /* Deep copy into result */
            result.words[result.count].len = cw.len;
            result.words[result.count].alphabet = cw.alphabet;
            result.words[result.count].symbols =
                (size_t *)malloc(cw.len * sizeof(size_t));
            if (result.words[result.count].symbols) {
                memcpy(result.words[result.count].symbols, cw.symbols,
                       cw.len * sizeof(size_t));
                result.count++;
            }
        }

        ld_free_codeword(&cw);
    }

    free(msg);
    free(candidate.symbols);
    return result;
}

/* ==================================================================
 *  L4 — Johnson Bound
 *  Theorem (Johnson 1962):
 *    For a code C ⊆ [q]^n with distance d,
 *    |B(x, e) ∩ C| ≤ q·n    for e ≤ n·(1 - sqrt{1 - d/n}).
 * ================================================================== */

size_t ld_johnson_bound(size_t n, size_t d, size_t q)
{
    if (n == 0 || d > n) return 0;

    double frac = (double)d / (double)n;
    double inner = 1.0 - frac;
    if (inner < 0.0) inner = 0.0;

    double sqrt_term = sqrt(inner);
    double e = (double)n * (1.0 - sqrt_term);

    /* Johnson bound list size: L_J ≤ q·n  for e ≤ e_J */
    size_t e_floor = (size_t)floor(e);

    /* Actually, the Johnson bound says: if e ≤ n·(1 - sqrt{1 - d/n}),
     * then list size is ≤ q·n.  We return the bound on list size. */
    return q * n;
}

double ld_johnson_radius(size_t n, size_t d)
{
    if (n == 0 || d > n) return 0.0;

    double frac = (double)d / (double)n;
    double inner = 1.0 - frac;
    if (inner < 0.0) inner = 0.0;

    return (double)n * (1.0 - sqrt(inner));
}

/* Combinatorial list-decoding capacity:
 * For an (n,k) code over alphabet q, the maximum error fraction δ
 * that can be list-decoded with polynomial list size is:
 *
 *   δ_cap = 1 - R   (where R = k/n is the rate)
 *
 * More precisely, for any ε > 0, there exist codes of rate
 * R = 1 - δ - ε that are list-decodable up to error fraction δ
 * with list size O(1/ε).
 *
 * Reference: Zyablov-Pinsker (1981), Elias (1991).
 */
double ld_list_decoding_capacity(double rate, double epsilon)
{
    /* Capacity: δ < 1 - R - ε */
    double delta_cap = 1.0 - rate - epsilon;
    return (delta_cap > 0.0) ? delta_cap : 0.0;
}

/* Griesmer bound for linear [n,k,d]_q codes:
 *
 *   n ≥ Σ_{i=0}^{k-1} ⌈ d / q^i ⌉
 *
 * This is a necessary condition for the existence of a linear code
 * with given parameters.  It's a strengthening of the Singleton bound.
 *
 * Reference: Griesmer (1960).
 */
size_t ld_griesmer_bound(size_t k, size_t d, size_t q)
{
    size_t n_sum = 0;
    size_t q_pow = 1;  /* q^0 = 1 */

    for (size_t i = 0; i < k; i++) {
        /* ⌈ d / q^i ⌉ */
        size_t term = (d + q_pow - 1) / q_pow;  /* ceiling division */
        n_sum += term;
        q_pow *= q;
    }

    return n_sum;
}

/* ==================================================================
 *  L2 — Helper Functions
 * ================================================================== */

bool ld_is_list_decodable(const CodeParams *cp, double radius, size_t L)
{
    if (!cp) return false;
    if (radius < 0.0 || radius > 1.0) return false;

    /* Compute Johnson bound list size for this radius */
    size_t e = (size_t)(radius * (double)cp->n);
    size_t johnson_L = ld_johnson_bound(cp->n, cp->d, cp->q);

    /* Check if list size L is sufficient */
    if (L >= johnson_L) {
        /* Well within Johnson bound — list-decodable */
        return true;
    }

    /* Check against capacity */
    double rate = (double)cp->k / (double)cp->n;
    double cap = 1.0 - rate;
    if (radius <= cap) {
        /* Within capacity — list-decodable with polynomial list size */
        return true;
    }

    return false;
}

void ld_free_codeword(Codeword *cw)
{
    if (!cw) return;
    free(cw->symbols);
    cw->symbols = NULL;
    cw->len = 0;
}

void ld_free_list(CodewordList *cl)
{
    if (!cl) return;
    for (size_t i = 0; i < cl->count; i++) {
        free(cl->words[i].symbols);
    }
    free(cl->words);
    cl->words = NULL;
    cl->count = 0;
}

Codeword ld_copy_codeword(const Codeword *src)
{
    Codeword cpy;
    cpy.len = 0;
    cpy.alphabet = 0;
    cpy.symbols = NULL;

    if (!src || !src->symbols || src->len == 0) return cpy;

    cpy.len = src->len;
    cpy.alphabet = src->alphabet;
    cpy.symbols = (size_t *)malloc(src->len * sizeof(size_t));
    if (cpy.symbols) {
        memcpy(cpy.symbols, src->symbols, src->len * sizeof(size_t));
    } else {
        cpy.len = 0;
    }

    return cpy;
}

void ld_print_codeword(const Codeword *cw)
{
    if (!cw || !cw->symbols) {
        printf("[empty codeword]\n");
        return;
    }

    printf("[");
    for (size_t i = 0; i < cw->len; i++) {
        printf("%zu", cw->symbols[i]);
        if (i + 1 < cw->len) printf(" ");
    }
    printf("] (len=%zu, q=%zu)\n", cw->len, cw->alphabet);
}

/* ==================================================================
 *  L5 — Brute-Force Codeword Enumeration (reference)
 * ================================================================== */

/**
 * Enumerate all codewords of a linear [n,k]_q code given
 * the generator matrix G (k×n).  Used as ground truth for
 * verifying list-decoding algorithms.
 *
 * @param G       Generator matrix (flat, k×n, row-major).
 * @param n       Block length.
 * @param k       Dimension.
 * @param q       Alphabet size.
 * @param cw_out  Output: array of q^k codewords (each length n).
 * @return Number of codewords (q^k).
 */
size_t ld_enumerate_codewords(const size_t *G, size_t n, size_t k, size_t q,
                               size_t **cw_out)
{
    /* Compute total number of codewords */
    size_t total = 1;
    for (size_t i = 0; i < k && total <= (size_t)1e5; i++) {
        total *= q;
    }
    if (total > (size_t)1e5) {
        *cw_out = NULL;
        return 0;  /* Too many codewords */
    }

    size_t *codewords = (size_t *)calloc(total * n, sizeof(size_t));
    if (!codewords) { *cw_out = NULL; return 0; }

    /* Enumerate all k-tuples over [q] */
    for (size_t idx = 0; idx < total; idx++) {
        /* Extract message digits (base q) */
        size_t *msg = (size_t *)calloc(k, sizeof(size_t));
        size_t tmp = idx;
        for (size_t i = 0; i < k; i++) {
            msg[i] = tmp % q;
            tmp /= q;
        }

        /* Compute codeword = msg × G  (over GF(q)) */
        for (size_t j = 0; j < n; j++) {
            size_t val = 0;
            for (size_t i = 0; i < k; i++) {
                val = (val + msg[i] * G[i * n + j]) % q;
            }
            codewords[idx * n + j] = val;
        }
        free(msg);
    }

    *cw_out = codewords;
    return total;
}
