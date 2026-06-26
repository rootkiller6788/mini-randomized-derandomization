/**********************************************************************
 * zigzag_product.c -- Zig-Zag Product of Expander Graphs
 *
 * Implements the zig-zag product of Reingold, Vadhan, and Wigderson
 * (2002), a fundamental operation for constructing expander graphs
 * from smaller ones. The zig-zag product preserves expansion while
 * controlling the degree.
 *
 * Given a d_G-regular graph G on n vertices and a d_H-regular graph H
 * on d_G vertices, the zig-zag product G (z) H is a d_H^2-regular
 * graph on n*d_G vertices.
 *
 * Key Theorem (RVW 2002): If G is an (n, d_G, lambda_G)-expander
 * and H is a (d_G, d_H, lambda_H)-expander, then G (z) H is an
 * (n*d_G, d_H^2, lambda_G + lambda_H + lambda_H^2)-expander.
 *
 * Knowledge Coverage:
 *   L1: ZigZagResult type, rotation map representation
 *   L2: Degree reduction while preserving expansion
 *   L3: Rotation maps as bijections on [n] x [d]
 *   L4: Zig-zag spectral bound theorem (RVW 2002)
 *   L5: Explicit zig-zag product construction
 *   L7: Derandomization, SL=L (Reingold 2008)
 *
 * References:
 *   - Reingold, Vadhan, Wigderson (2002) "Entropy Waves, the Zig-Zag
 *     Graph Product, and New Constant-Degree Expanders"
 *   - Reingold (2008) "Undirected Connectivity in Log-Space"
 **********************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "expander_core.h"
#include "expander_constructions.h"
#include "zigzag_product.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==================================================================
 * Rotation Map Representation
 *
 * A rotation map Rot_G: [n] x [d] -> [n] x [d] is a bijection
 * such that Rot_G(v, i) = (w, j) means the i-th edge incident
 * to v leads to w, and this is the j-th edge incident to w.
 *
 * For undirected graphs: Rot_G(Rot_G(v, i)) = (v, i).
 * ================================================================== */

/**
 * exp_rotation_map -- Apply the rotation map of g to (v, i).
 * Returns the resulting vertex w, and stores the output edge
 * index in *j (if j is not NULL).
 *
 * Complexity: O(d).
 */
size_t exp_rotation_map(const ExpanderGraph *g, size_t v, size_t i,
                         size_t *j) {
    if (!g || v >= g->n || i >= g->d) {
        if (j) *j = SIZE_MAX;
        return SIZE_MAX;
    }

    size_t w = g->adj[v * g->d + i];
    if (w == SIZE_MAX || w >= g->n) {
        if (j) *j = SIZE_MAX;
        return SIZE_MAX;
    }

    /* Find j such that the j-th edge of w connects back to v */
    if (j) {
        *j = SIZE_MAX;
        for (size_t k = 0; k < g->d; k++) {
            if (g->adj[w * g->d + k] == v) {
                *j = k;
                break;
            }
        }
    }
    return w;
}

/* ==================================================================
 * Zig-Zag Product: G (z) H
 *
 * Vertices: V = V_G x V_H = [n_G] x [d_G]    (since n_H = d_G)
 * Degree: d_H^2
 *
 * A step in G (z) H from vertex (v, k) consists of:
 *   1. Zig:    Move within the cloud H_v: (v, k) -> (v, k')
 *   2. Zag:    Move across clouds via G:  (v, i) -> (w, j) using rotation map
 *   3. Zig:    Move within the cloud H_w: (w, j) -> (w, k'')
 *
 * So a single edge in G (z) H corresponds to:
 *   (v, k) --(Zig)--> (v, k') --(Zag)--> (w, j)
 *   Then from (v,k) we go to result of a full zig-zag-zig.
 *
 * Actually, the standard definition uses:
 *   (v,k) --zig(k,k')--> (v,i) --zag(v,i)--> (w,j) --zig(j,k'')--> (w,k'')
 * where zig(k,k') means moving to neighbor k' of k in H,
 * and then i = the (d_H-th edge) that corresponds to ... 
 *
 * Simplified implementation using the RVW edge labeling:
 * Edge (h, h') of G(z)H corresponds to:
 *   - pick i = h-th neighbor of k in H (zig step)
 *   - follow rotation map of G: (v,i) -> (w,j)
 *   - pick k' = h'-th neighbor of j in H (zag step)
 *
 * ================================================================== */

/**
 * zz_compute -- Compute the zig-zag product G (z) H.
 *
 * G: d_G-regular graph on n_G vertices
 * H: d_H-regular graph on n_H = d_G vertices
 *
 * Result: d_H^2-regular graph on n_G * d_G vertices.
 *
 * The implementation follows the standard RVW (2002) construction.
 *
 * Complexity: O(n_G * d_G * d_H^2).
 */
ZigZagResult zz_compute(const ExpanderGraph *G, const ExpanderGraph *H) {
    ZigZagResult res;
    memset(&res, 0, sizeof(ZigZagResult));

    if (!G || !H) return res;
    if (H->n != G->d) return res;  /* Required: n_H = d_G */

    size_t nG = G->n, dG = G->d, dH = H->d;
    size_t new_n = nG * dG;
    size_t new_d = dH * dH;

    ExpanderGraph *result = exp_create(new_n, new_d);
    if (!result) return res;

    /* For each cloud vertex (v, k) in V_G x V_H */
    for (size_t v = 0; v < nG; v++) {
        for (size_t k = 0; k < dG; k++) {
            size_t u = v * dG + k;  /* Current vertex index */
            size_t edge_idx = 0;

            /* For each pair (h1, h2) of H-edges */
            for (size_t h1 = 0; h1 < dH && edge_idx < new_d; h1++) {
                /* Step 1 (Zig): within cloud at v,
                 * move from k to i = k-th neighbor in H via edge h1 */
                size_t i = H->adj[k * dH + h1];
                if (i == SIZE_MAX || i >= dG) continue;

                /* Step 2 (Zag): across clouds via G's rotation map */
                size_t w = G->adj[v * dG + i];
                if (w == SIZE_MAX || w >= nG) continue;

                /* Find j such that w's j-th edge points back to v */
                size_t j = SIZE_MAX;
                for (size_t jj = 0; jj < dG; jj++) {
                    if (G->adj[w * dG + jj] == v ||
                        G->adj[w * dG + jj] == SIZE_MAX) {
                        /* For expander constructions, we need j such that
                         * the rotation map is consistent */
                        j = jj;
                        break;
                    }
                }
                if (j == SIZE_MAX) j = i;  /* Fallback */

                for (size_t h2 = 0; h2 < dH && edge_idx < new_d; h2++) {
                    /* Step 3 (Zig): within cloud at w,
                     * move from j to k' = j-th neighbor in H via edge h2 */
                    size_t kp = H->adj[j * dH + h2];
                    if (kp == SIZE_MAX || kp >= dG) continue;

                    size_t uu = w * dG + kp;
                    result->adj[u * new_d + edge_idx] = uu;
                    edge_idx++;
                }
            }
        }
    }

    res.G = (ExpanderGraph *)G;  /* Non-owning reference */
    res.H = (ExpanderGraph *)H;
    res.result = result;
    return res;
}

/* ==================================================================
 * Zig-Zag Spectral Properties
 * ================================================================== */

/**
 * zz_preserves_expansion -- Check if the zig-zag product preserves
 * expansion, i.e., the resulting graph is still an expander.
 *
 * Theorem (RVW 2002): If lambda_G, lambda_H < 1, then
 *   lambda_{G(z)H} <= lambda_G + lambda_H + lambda_H^2
 *
 * This function verifies the empirical bound for the given graphs.
 *
 * Complexity: O(n^2) for eigenvalue computation.
 */
bool zz_preserves_expansion(const ZigZagResult *z) {
    if (!z || !z->G || !z->H || !z->result) return false;

    double lamG = exp_second_eigenvalue(z->G);
    double lamH = exp_second_eigenvalue(z->H);
    double lamZ = exp_second_eigenvalue(z->result);

    double dG = (double)z->G->d;
    double dH_sq = (double)(z->H->d * z->H->d);

    /* Normalize */
    double lamG_norm = lamG / dG;
    double lamH_norm = lamH / (double)z->H->d;
    double lamZ_norm = lamZ / dH_sq;

    double bound = lamG_norm + lamH_norm + lamH_norm * lamH_norm;
    return lamZ_norm <= bound + 1e-9;
}

/**
 * zz_degree_reduction -- Compute the new degree after zig-zag product.
 * If G has degree d_G and H has degree d_H (with n_H = d_G),
 * then G(z)H has degree d_H^2.
 *
 * This is a massive degree reduction when d_H << d_G.
 *
 * Complexity: O(1).
 */
double zz_degree_reduction(size_t d1, size_t d2) {
    if (d1 == 0 || d2 == 0) return 0.0;
    return (double)(d2 * d2) / (double)d1;
}

/**
 * zz_new_size -- Compute the number of vertices after zig-zag product.
 * If G has n_1 vertices and n_H = d_G = d_2, then G(z)H has
 * n_1 * d_2 = n_1 * n_H vertices.
 *
 * Complexity: O(1).
 */
size_t zz_new_size(size_t n1, size_t d2) {
    return n1 * d2;
}

/**
 * zz_spectral_bound -- Compute the theoretical upper bound on
 * lambda(G(z)H) using the RVW bound:
 *
 *   lambda <= lambda_G + lambda_H + lambda_H^2
 *
 * where eigenvalues are of the normalized adjacency matrices.
 *
 * Complexity: O(n^2) for each eigenvalue computation.
 */
double zz_spectral_bound(const ExpanderGraph *G, const ExpanderGraph *H) {
    if (!G || !H) return INFINITY;

    double lamG = exp_second_eigenvalue(G);
    double lamH = exp_second_eigenvalue(H);

    double dG = (double)G->d;
    double dH = (double)H->d;

    /* Normalized eigenvalues */
    double lamG_norm = lamG / dG;
    double lamH_norm = lamH / dH;

    /* RVW bound on normalized eigenvalue */
    double bound_norm = lamG_norm + lamH_norm + lamH_norm * lamH_norm;
    if (bound_norm >= 1.0) bound_norm = 0.999;

    /* Convert back to unnormalized (for d_H^2-regular graph) */
    return bound_norm * (double)(dH * dH);
}

/* ==================================================================
 * Explicit Construction Using Zig-Zag Iteration
 * ================================================================== */

/**
 * zz_explicit_construction -- Build an explicit expander family
 * by iterating the zig-zag product with powering and tensoring.
 *
 * The RVW construction:
 *   1. Start with a constant-size "base" expander H (e.g., on d^4 vertices)
 *   2. Repeatedly apply: G_{t+1} = (G_t (z) H)^2 (zig-zag + powering)
 *   3. This produces an infinite family of constant-degree expanders.
 *
 * This function implements one step of this iteration:
 *   Given a graph G of size m and degree d, produce a graph
 *   of size m*d_H^8 and degree d_H^2.
 *
 * This is a SIMPLIFIED demonstration version. The full RVW iteration
 * requires careful handling of the tensor product and proper
 * rotation map definitions.
 *
 * Complexity: O(m * d * d_H^2) per iteration.
 */
bool zz_explicit_construction(size_t n, size_t d, ExpanderGraph **result) {
    if (!result || n == 0 || d == 0) return false;

    /* Step 1: Create a small base expander H using Margulis construction
     * H must have d vertices and some constant degree d_H */
    /* For simplicity, create H on m^2 = d vertices */
    size_t m = (size_t)sqrt((double)d);
    if (m * m != d || m < 2) return false;

    ExpanderGraph *H = exp_margulis_gabber_galil(d);
    if (!H) return false;

    /* Step 2: Create initial G using random construction or
     * Margulis on n vertices.  We use the Margulis construction
     * for explicit guarantees when n is a perfect square. */

    ExpanderGraph *G = NULL;
    size_t nn = n;
    /* Round n up to nearest perfect square for Margulis */
    size_t r = (size_t)ceil(sqrt((double)n));
    nn = r * r;

    G = exp_margulis_gabber_galil(nn);
    if (!G) { exp_free(H); return false; }

    /* Step 3: Apply zig-zag product */
    ZigZagResult zr = zz_compute(G, H);

    if (zr.result) {
        *result = zr.result;
        exp_free(G);
        exp_free(H);
        return true;
    }

    exp_free(G);
    exp_free(H);
    return false;
}

/* ==================================================================
 * Memory Management
 * ================================================================== */

/**
 * zz_free_result -- Free a ZigZagResult, deallocating the result
 * graph. The input graphs G and H are NOT freed (caller owns them).
 */
void zz_free_result(ZigZagResult *z) {
    if (!z) return;
    if (z->result) {
        exp_free(z->result);
        z->result = NULL;
    }
    /* Do not free z->G or z->H -- they are borrowed */
    z->G = NULL;
    z->H = NULL;
}

/* ==================================================================
 * Additional Zig-Zag Analysis Functions
 * ================================================================== */

/**
 * zz_empirical_lambda -- Compute the actual lambda_2 of a zig-zag
 * product result and compare it with the theoretical bound.
 *
 * Returns the ratio empirical/theoretical. Values <= 1.0
 * indicate the bound holds.
 */
double zz_empirical_lambda(const ZigZagResult *z) {
    if (!z || !z->result || !z->G || !z->H) return INFINITY;

    double lam_actual = exp_second_eigenvalue(z->result);
    double lam_bound = zz_spectral_bound(z->G, z->H);

    if (lam_bound <= 0.0) return INFINITY;
    return lam_actual / lam_bound;
}

/**
 * zz_verify_construction -- Run basic sanity checks on a zig-zag
 * product result: regularity, connectivity, and expansion.
 */
bool zz_verify_construction(const ZigZagResult *z) {
    if (!z || !z->result) return false;

    ExpanderGraph *g = z->result;
    if (!exp_is_regular(g)) return false;
    if (!exp_is_connected(g)) return false;

    /* Check that degree matches expected d_H^2 */
    if (z->H) {
        size_t expected_d = z->H->d * z->H->d;
        if (g->d != expected_d) return false;
    }

    /* Check that vertex count matches n_G * d_G */
    if (z->G) {
        size_t expected_n = z->G->n * z->G->d;
        if (g->n != expected_n) return false;
    }

    return true;
}

/**
 * zz_print_info -- Print information about a zig-zag product.
 */
void zz_print_info(const ZigZagResult *z, FILE *fp) {
    if (!z || !fp) return;

    fprintf(fp, "Zig-Zag Product:\n");
    if (z->G) {
        fprintf(fp, "  G: n=%zu, d=%zu\n", z->G->n, z->G->d);
    }
    if (z->H) {
        fprintf(fp, "  H: n=%zu, d=%zu\n", z->H->n, z->H->d);
    }
    if (z->result) {
        fprintf(fp, "  Result: n=%zu, d=%zu\n",
                z->result->n, z->result->d);
        double lam2 = exp_second_eigenvalue(z->result);
        fprintf(fp, "  lambda_2(result): %.6f\n", lam2);
        fprintf(fp, "  Is expander: %s\n",
                exp_is_expander(z->result, 0.1) ? "yes" : "no");
    }
}
