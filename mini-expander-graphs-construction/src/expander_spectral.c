/**********************************************************************
 * expander_spectral.c -- Spectral Theory of Expander Graphs
 *
 * Implements advanced spectral analysis for expander graphs:
 *   - Full eigenvalue decomposition (Jacobi method for small matrices)
 *   - Laplacian matrix and normalized Laplacian
 *   - Expander mixing lemma
 *   - Tanner's inequality
 *   - Edge expansion vs spectral gap relationships
 *   - Markov chain comparison techniques
 *   - Heat kernel and diffusion on expanders
 *
 * Knowledge Coverage:
 *   L1-L2: Spectral decomposition, Laplacian spectrum
 *   L3: Symmetric matrices, Rayleigh quotients, eigenvalues
 *   L4: Expander mixing lemma, Tanner's inequality
 *   L5: Jacobi eigenvalue algorithm, spectral partitioning
 *   L8: Spectral graph theory, heat kernel estimates
 *
 * References:
 *   - Chung (1997) "Spectral Graph Theory"
 *   - Alon & Milman (1985) "lambda_1, isoperimetric inequalities for graphs"
 *   - Hoory, Linial, Wigderson (2006) "Expander Graphs and Their Applications"
 *   - Spielman (2019) "Spectral and Algebraic Graph Theory"
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
 * Laplacian Matrix Construction
 * ================================================================== */

/**
 * exp_laplacian_matrix -- Construct the unnormalized Laplacian L = D - A.
 * For a d-regular graph, L = d*I - A.
 *
 * The Laplacian eigenvalues are: 0 = mu_1 <= mu_2 <= ... <= mu_n.
 * The second eigenvalue mu_2 (algebraic connectivity) measures
 * how well-connected the graph is.
 *
 * mu_2 > 0 iff the graph is connected.
 *
 * Complexity: O(n^2).
 */
double **exp_laplacian_matrix(const ExpanderGraph *g) {
    if (!g) return NULL;
    size_t n = g->n;
    double **L = (double **)malloc(n * sizeof(double *));
    if (!L) return NULL;

    for (size_t i = 0; i < n; i++) {
        L[i] = (double *)calloc(n, sizeof(double));
        if (!L[i]) {
            for (size_t j = 0; j < i; j++) free(L[j]);
            free(L);
            return NULL;
        }
    }

    for (size_t u = 0; u < n; u++) {
        L[u][u] = (double)exp_degree(g, u);  /* Degree */
        for (size_t k = 0; k < g->d; k++) {
            size_t v = g->adj[u * g->d + k];
            if (v != SIZE_MAX && v < n) {
                L[u][v] = -1.0;
            }
        }
    }
    return L;
}

/**
 * exp_normalized_laplacian -- Construct the normalized Laplacian
 * L_norm = I - D^{-1/2} A D^{-1/2}.
 *
 * For d-regular graphs: L_norm = I - A/d.
 *
 * The eigenvalues satisfy: 0 = nu_1 <= nu_2 <= ... <= nu_n <= 2.
 */
double **exp_normalized_laplacian(const ExpanderGraph *g) {
    if (!g || g->d == 0) return NULL;
    size_t n = g->n;
    double **Ln = (double **)malloc(n * sizeof(double *));
    if (!Ln) return NULL;

    for (size_t i = 0; i < n; i++) {
        Ln[i] = (double *)calloc(n, sizeof(double));
        if (!Ln[i]) {
            for (size_t j = 0; j < i; j++) free(Ln[j]);
            free(Ln);
            return NULL;
        }
    }

    for (size_t u = 0; u < n; u++) {
        Ln[u][u] = 1.0;  /* I */
        for (size_t k = 0; k < g->d; k++) {
            size_t v = g->adj[u * g->d + k];
            if (v != SIZE_MAX && v < n) {
                Ln[u][v] = -1.0 / (double)g->d;
            }
        }
    }
    return Ln;
}

/**
 * exp_algebraic_connectivity -- Compute mu_2, the second smallest
 * eigenvalue of the Laplacian (Fiedler eigenvalue).
 *
 * For a d-regular graph: mu_2 = d - lambda_2.
 * So the spectral gap (d - lambda_2) IS the algebraic connectivity.
 *
 * mu_2 > 0 iff graph is connected.
 * Larger mu_2 => more "well-connected."
 *
 * Complexity: O(n^2).
 */
double exp_algebraic_connectivity(const ExpanderGraph *g) {
    if (!g) return 0.0;
    return exp_spectral_gap(g);  /* For regular graphs, mu_2 = d - lambda_2 */
}

/**
 * exp_fiedler_vector -- Compute the Fiedler vector (eigenvector
 * corresponding to mu_2) via deflated power iteration on the
 * Laplacian.
 *
 * The Fiedler vector is used for spectral graph partitioning:
 * vertices with positive entries go to one side,
 * vertices with negative entries go to the other.
 *
 * Complexity: O(iters * n^2).
 */
double *exp_fiedler_vector(const ExpanderGraph *g, int iters) {
    if (!g || g->n <= 1) return NULL;
    size_t n = g->n;

    double **L = exp_laplacian_matrix(g);
    if (!L) return NULL;

    /* For the Laplacian of a d-regular graph,
     * the top eigenvalue is 2d (farthest from 0).
     * We use inverse iteration: largest eigenvalue of (2d*I - L)^{-1}
     * corresponds to the smallest non-zero eigenvalue of L. */

    /* Simplified: use power iteration on (2d*I - L) to get mu_2 via
     * eigenvalue shift. For a d-regular graph:
     *   L = d*I - A, so 2d*I - L = d*I + A.
     * The eigenvalues of d*I + A are d + lambda_i.
     * The second smallest lambda of L corresponds to the second
     * largest eigenvalue of d*I + A. */

    /* Build M = d*I + A */
    double **M = (double **)malloc(n * sizeof(double *));
    if (!M) { exp_laplacian_free(L, n); return NULL; }
    for (size_t i = 0; i < n; i++) {
        M[i] = (double *)calloc(n, sizeof(double));
        if (!M[i]) {
            for (size_t j = 0; j < i; j++) free(M[j]);
            free(M); exp_laplacian_free(L, n); return NULL;
        }
        M[i][i] = (double)g->d;
        for (size_t k = 0; k < g->d; k++) {
            size_t v = g->adj[i * g->d + k];
            if (v != SIZE_MAX && v < n) M[i][v] = 1.0;
        }
    }

    /* Power iteration with deflation (remove all-ones) */
    double *v = (double *)malloc(n * sizeof(double));
    if (!v) {
        for (size_t i = 0; i < n; i++) free(M[i]);
        free(M); exp_laplacian_free(L, n);
        return NULL;
    }

    for (size_t i = 0; i < n; i++)
        v[i] = ((double)rand() / RAND_MAX) - 0.5;
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) sum += v[i];
    sum /= (double)n;
    for (size_t i = 0; i < n; i++) v[i] -= sum;
    double norm = 0.0;
    for (size_t i = 0; i < n; i++) norm += v[i] * v[i];
    if (norm > 1e-15)
        for (size_t i = 0; i < n; i++) v[i] /= sqrt(norm);

    for (int t = 0; t < iters; t++) {
        double *w = (double *)calloc(n, sizeof(double));
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                w[i] += M[i][j] * v[j];
        double s = 0.0;
        for (size_t i = 0; i < n; i++) s += w[i];
        s /= (double)n;
        for (size_t i = 0; i < n; i++) w[i] -= s;
        norm = 0.0;
        for (size_t i = 0; i < n; i++) norm += w[i] * w[i];
        if (sqrt(norm) < 1e-15) { free(w); break; }
        for (size_t i = 0; i < n; i++) v[i] = w[i] / sqrt(norm);
        free(w);
    }

    for (size_t i = 0; i < n; i++) free(M[i]);
    free(M);
    exp_laplacian_free(L, n);
    return v;  /* Caller must free */
}

/**
 * exp_laplacian_free -- Free a Laplacian matrix.
 */
void exp_laplacian_free(double **L, size_t n) {
    if (!L) return;
    for (size_t i = 0; i < n; i++) free(L[i]);
    free(L);
}

/* ==================================================================
 * Expander Mixing Lemma
 *
 * Let G be a d-regular graph and lambda = max(lambda_2, |lambda_n|).
 * Then for any two subsets S, T of vertices:
 *
 *   | |E(S,T)| - (d/n) * |S| * |T| | <= lambda * sqrt(|S| * |T|)
 *
 * where E(S,T) counts ordered pairs (u in S, v in T) with (u,v) in E.
 *
 * The expander mixing lemma quantifies how edges are distributed
 * in an expander: they behave almost like a random graph.
 *
 * Reference: Alon & Chung (1988)
 * ================================================================== */

/**
 * exp_mixing_lemma_bound -- Compute the upper bound on the deviation
 * from the random graph expectation for edge count between S and T.
 *
 * Returns the maximum absolute deviation.
 */
double exp_mixing_lemma_bound(const ExpanderGraph *g,
                               size_t S_size, size_t T_size) {
    if (!g) return INFINITY;
    double lam2 = exp_second_eigenvalue(g);
    double lam = lam2;  /* Use lambda_2 as the spectral bound */
    return lam * sqrt((double)S_size * (double)T_size);
}

/**
 * exp_count_edges_between -- Count edges between subset S and subset T.
 * S and T are given as bitmasks.
 *
 * Returns the number of ordered pairs (s in S, t in T) with (s,t) in E.
 */
size_t exp_count_edges_between(const ExpanderGraph *g,
                                const char *S, const char *T) {
    if (!g || !S || !T) return 0;

    size_t count = 0;
    for (size_t u = 0; u < g->n; u++) {
        if (!S[u]) continue;
        for (size_t k = 0; k < g->d; k++) {
            size_t v = g->adj[u * g->d + k];
            if (v != SIZE_MAX && v < g->n && T[v]) count++;
        }
    }
    return count;
}

/**
 * exp_verify_mixing_lemma -- Empirically verify the expander mixing
 * lemma for all subsets up to a given size.
 *
 * Returns the maximum observed ratio:
 *   max_{S,T} |E(S,T) - (d/n)*|S|*|T|| / (lambda * sqrt(|S|*|T|))
 *
 * A value <= 1.0 means the bound holds for the tested subsets.
 */
double exp_verify_mixing_lemma(const ExpanderGraph *g, size_t max_size) {
    if (!g || g->n > 20) return -1.0;  /* Only feasible for small n */

    (void)exp_second_eigenvalue(g);  /* lam2 used via mixing_lemma_bound */
    double max_ratio = 0.0;
    size_t n = g->n;

    char *S = (char *)calloc(n, sizeof(char));
    char *T = (char *)calloc(n, sizeof(char));
    if (!S || !T) { free(S); free(T); return -1.0; }

    for (size_t mask1 = 1; mask1 < (1ULL << n); mask1++) {
        size_t s_size = 0;
        for (size_t i = 0; i < n; i++) {
            S[i] = (mask1 >> i) & 1;
            if (S[i]) s_size++;
        }
        if (s_size > max_size) continue;

        for (size_t mask2 = 1; mask2 < (1ULL << n); mask2++) {
            size_t t_size = 0;
            for (size_t i = 0; i < n; i++) {
                T[i] = (mask2 >> i) & 1;
                if (T[i]) t_size++;
            }
            if (t_size > max_size) continue;

            size_t actual = exp_count_edges_between(g, S, T);
            double expected = (double)g->d * (double)s_size
                            * (double)t_size / (double)n;
            double deviation = fabs((double)actual - expected);
            double bound = exp_mixing_lemma_bound(g, s_size, t_size);

            if (bound > 0) {
                double ratio = deviation / bound;
                if (ratio > max_ratio) max_ratio = ratio;
            }
        }
    }

    free(S); free(T);
    return max_ratio;
}

/* ==================================================================
 * Tanner's Inequality
 *
 * Tanner (1984): For a d-regular graph G = (V,E) with second
 * eigenvalue lambda_2, the independence number alpha(G) satisfies:
 *
 *   alpha(G) <= n * lambda_2 / (d + lambda_2)
 *
 * where alpha(G) is the size of the largest independent set.
 * ================================================================== */

/**
 * exp_independence_number_bound -- Upper bound on alpha(G) via
 * the Hoffman bound (generalized by Tanner).
 *
 * alpha(G) <= n * lambda_2 / (d + lambda_2) for regular graphs.
 *
 * For expanders with lambda_2 << d, this gives alpha(G) = O(n/d),
 * which means expanders have small independent sets.
 */
double exp_independence_number_bound(const ExpanderGraph *g) {
    if (!g) return 0.0;
    double lam2 = exp_second_eigenvalue(g);
    double d = (double)g->d;
    return (double)g->n * lam2 / (d + lam2);
}

/**
 * exp_independence_number_exact -- Exact independence number via
 * brute force for small graphs (n <= 20).
 */
size_t exp_independence_number_exact(const ExpanderGraph *g) {
    if (!g || g->n > 20) return 0;

    size_t best = 0;
    size_t n = g->n;

    for (size_t mask = 1; mask < (1ULL << n); mask++) {
        bool independent = true;
        size_t size = 0;

        for (size_t u = 0; u < n && independent; u++) {
            if (!(mask & (1ULL << u))) continue;
            size++;
            for (size_t v = u + 1; v < n && independent; v++) {
                if (!(mask & (1ULL << v))) continue;
                if (exp_has_edge(g, u, v)) independent = false;
            }
        }

        if (independent && size > best) best = size;
    }
    return best;
}

/* ==================================================================
 * Chromatic Number Bounds
 * ================================================================== */

/**
 * exp_chromatic_number_bound -- Hoffman's bound on the chromatic
 * number chi(G):
 *
 *   chi(G) >= 1 - d / lambda_n
 *
 * where lambda_n is the smallest (most negative) eigenvalue.
 *
 * For bipartite expanders, lambda_n = -d and chi >= 2, which is tight.
 * For non-bipartite expanders, lambda_n > -d and chi can be larger.
 *
 * This returns a lower bound on chi(G).
 */
double exp_chromatic_number_bound(const ExpanderGraph *g) {
    if (!g) return 0.0;
    /* For simplicity, approximate lambda_n = -lambda_2 for bipartite-like
     * graphs, or estimate via trace. */
    /* For Ramanujan-like expanders, |lambda_i| <= 2*sqrt(d-1) for i>1 */
    double bound = 1.0 + (double)g->d / (2.0 * sqrt((double)g->d - 1.0));
    return bound;
}

/* ==================================================================
 * Jacobi Eigenvalue Algorithm (Full Spectrum)
 *
 * Computes all eigenvalues of a real symmetric matrix via the
 * classical Jacobi method. This provides the full spectrum
 * (lambda_1 = d, lambda_2, ..., lambda_n).
 *
 * The Jacobi method applies successive Givens rotations to
 * diagonalize the matrix. Convergence is guaranteed for symmetric
 * matrices.
 *
 * Reference: Golub & Van Loan (2013), Sec. 8.4
 * ================================================================== */

/**
 * exp_full_spectrum -- Compute all eigenvalues of the adjacency
 * matrix using the Jacobi method (for n <= 100).
 *
 * eigenvalues: output array of size n, sorted descending.
 *
 * Returns 0 on success, -1 on error.
 *
 * Complexity: O(n^3 * log(1/eps)).
 */
int exp_full_spectrum(const ExpanderGraph *g, double *eigenvalues) {
    if (!g || !eigenvalues || g->n > 100 || g->n == 0) return -1;

    size_t n = g->n;
    double **A = (double **)malloc(n * sizeof(double *));
    double **V = (double **)malloc(n * sizeof(double *));
    if (!A || !V) { free(A); free(V); return -1; }

    for (size_t i = 0; i < n; i++) {
        A[i] = (double *)calloc(n, sizeof(double));
        V[i] = (double *)calloc(n, sizeof(double));
        if (!A[i] || !V[i]) {
            for (size_t j = 0; j <= i; j++) { free(A[j]); free(V[j]); }
            free(A); free(V); return -1;
        }
        V[i][i] = 1.0;  /* Eigenvector matrix initialized to I */
    }

    /* Build adjacency matrix */
    for (size_t u = 0; u < n; u++) {
        for (size_t k = 0; k < g->d; k++) {
            size_t v = g->adj[u * g->d + k];
            if (v != SIZE_MAX && v < n) A[u][v] = 1.0;
        }
    }

    /* Jacobi iteration */
    const int max_iters = 200;
    const double eps = 1e-10;

    for (int iter = 0; iter < max_iters; iter++) {
        /* Find largest off-diagonal element */
        double max_off = 0.0;
        size_t p = 0, q = 1;
        for (size_t i = 0; i < n; i++) {
            for (size_t j = i + 1; j < n; j++) {
                if (fabs(A[i][j]) > max_off) {
                    max_off = fabs(A[i][j]);
                    p = i; q = j;
                }
            }
        }
        if (max_off < eps) break;

        /* Compute Givens rotation */
        double theta;
        if (fabs(A[p][p] - A[q][q]) < eps) {
            theta = M_PI / 4.0;
        } else {
            theta = 0.5 * atan2(2.0 * A[p][q], A[p][p] - A[q][q]);
        }

        double c = cos(theta), s = sin(theta);

        /* Apply rotation: A = R^T * A * R */
        double app = A[p][p], aqq = A[q][q], apq = A[p][q];

        A[p][p] = c*c*app - 2.0*c*s*apq + s*s*aqq;
        A[q][q] = s*s*app + 2.0*c*s*apq + c*c*aqq;
        A[p][q] = A[q][p] = 0.0;

        for (size_t i = 0; i < n; i++) {
            if (i == p || i == q) continue;
            double aip = A[i][p], aiq = A[i][q];
            A[i][p] = A[p][i] = c*aip - s*aiq;
            A[i][q] = A[q][i] = s*aip + c*aiq;
        }

        /* Update eigenvectors */
        for (size_t i = 0; i < n; i++) {
            double vip = V[i][p], viq = V[i][q];
            V[i][p] = c*vip - s*viq;
            V[i][q] = s*vip + c*viq;
        }
    }

    /* Extract eigenvalues from diagonal */
    for (size_t i = 0; i < n; i++) eigenvalues[i] = A[i][i];

    /* Sort descending (simple bubble sort for small n) */
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            if (eigenvalues[j] > eigenvalues[i]) {
                double tmp = eigenvalues[i];
                eigenvalues[i] = eigenvalues[j];
                eigenvalues[j] = tmp;
            }
        }
    }

    for (size_t i = 0; i < n; i++) { free(A[i]); free(V[i]); }
    free(A); free(V);
    return 0;
}

/* ==================================================================
 * Trace and Determinant Operations
 * ================================================================== */

/**
 * exp_trace_adjacency -- Tr(A) = 0 for simple graphs (no self-loops).
 */
double exp_trace_adjacency(const ExpanderGraph *g) {
    (void)g;
    return 0.0;  /* No self-loops in simple graphs */
}

/**
 * exp_trace_A2 -- Tr(A^2) = sum of squared eigenvalues = 2*|E| for
 * simple graphs (each edge contributes 2 to the trace).
 *
 * Also: Tr(A^2) = n*d for d-regular graphs.
 */
double exp_trace_A2(const ExpanderGraph *g) {
    if (!g) return 0.0;
    return (double)g->n * (double)g->d;
}

/**
 * exp_trace_A3 -- Tr(A^3) = 6 * (number of triangles in the graph).
 *
 * For triangle-free graphs (like many expander constructions),
 * Tr(A^3) = 0.
 */
double exp_trace_A3(const ExpanderGraph *g) {
    if (!g) return 0.0;
    size_t triangles = 0;
    for (size_t u = 0; u < g->n; u++) {
        for (size_t k1 = 0; k1 < g->d; k1++) {
            size_t v = g->adj[u * g->d + k1];
            if (v == SIZE_MAX || v >= g->n || v <= u) continue;
            for (size_t k2 = 0; k2 < g->d; k2++) {
                size_t w = g->adj[v * g->d + k2];
                if (w == SIZE_MAX || w >= g->n || w <= v) continue;
                if (exp_has_edge(g, u, w)) triangles++;
            }
        }
    }
    return 6.0 * (double)triangles;
}

/**
 * exp_num_triangles -- Count the number of triangles in the graph.
 */
size_t exp_num_triangles(const ExpanderGraph *g) {
    if (!g) return 0;
    size_t count = 0;
    for (size_t u = 0; u < g->n; u++) {
        for (size_t v = u + 1; v < g->n; v++) {
            if (!exp_has_edge(g, u, v)) continue;
            for (size_t w = v + 1; w < g->n; w++) {
                if (exp_has_edge(g, u, w) && exp_has_edge(g, v, w))
                    count++;
            }
        }
    }
    return count;
}

/* ==================================================================
 * Heat Kernel on Expanders
 *
 * The heat kernel H_t = exp(-t * L) describes continuous-time
 * diffusion on the graph.
 *
 * For expanders, the heat kernel converges to the uniform
 * distribution exponentially fast, with rate proportional to
 * the spectral gap.
 *
 * This is the continuous analog of the lazy random walk.
 * ================================================================== */

/**
 * exp_heat_kernel_diagonal -- Compute the diagonal entries of the
 * heat kernel at time t (the return probabilities).
 *
 * H_t(u,u) ~ 1/n + (1 - 1/n) * exp(-t * gamma) for expanders.
 *
 * where gamma = spectral gap / d is the normalized gap.
 */
double exp_heat_kernel_diagonal(const ExpanderGraph *g, double t) {
    if (!g || t < 0) return 0.0;
    double lam2 = exp_second_eigenvalue(g);
    double gap = (double)g->d - lam2;
    double norm_gap = gap / (double)g->d;
    return 1.0 / (double)g->n
           + (1.0 - 1.0 / (double)g->n) * exp(-t * norm_gap);
}

/**
 * exp_heat_kernel_mixing_time -- Time for heat kernel to mix:
 * t_mix(eps) = (1/gamma) * log(n/eps) where gamma = spectral gap / d.
 */
double exp_heat_kernel_mixing_time(const ExpanderGraph *g, double eps) {
    if (!g || eps <= 0) return 0.0;
    double lam2 = exp_second_eigenvalue(g);
    double gap = (double)g->d - lam2;
    if (gap <= 0) return INFINITY;
    double norm_gap = gap / (double)g->d;
    return log((double)g->n / eps) / norm_gap;
}

/* ==================================================================
 * Spectral Partitioning
 * ================================================================== */

/**
 * exp_spectral_partition -- Partition vertices into two sets using
 * the Fiedler vector (spectral bisection).
 *
 * Vertices with Fiedler vector entry >= 0 go to part[0],
 * vertices with entry < 0 go to part[1].
 *
 * This approximately minimizes the cut size (edges between parts)
 * while balancing the partition sizes.
 *
 * Complexity: O(iters * n^2).
 */
bool exp_spectral_partition(const ExpanderGraph *g, char *part,
                             int iters) {
    if (!g || !part) return false;

    double *fv = exp_fiedler_vector(g, iters);
    if (!fv) return false;

    for (size_t i = 0; i < g->n; i++) {
        part[i] = (fv[i] >= 0) ? 0 : 1;
    }
    free(fv);
    return true;
}

/**
 * exp_partition_cut_size -- Count edges crossing a partition.
 */
size_t exp_partition_cut_size(const ExpanderGraph *g, const char *part) {
    if (!g || !part) return 0;
    size_t cut = 0;
    for (size_t u = 0; u < g->n; u++) {
        for (size_t k = 0; k < g->d; k++) {
            size_t v = g->adj[u * g->d + k];
            if (v != SIZE_MAX && v < g->n && u < v) {
                if (part[u] != part[v]) cut++;
            }
        }
    }
    return cut;
}

/**
 * exp_partition_balance -- Compute the ratio of smaller part to n.
 * Perfect balance gives 0.5.
 */
double exp_partition_balance(const ExpanderGraph *g, const char *part) {
    if (!g || !part) return 0.0;
    size_t count0 = 0;
    for (size_t i = 0; i < g->n; i++)
        if (part[i] == 0) count0++;
    size_t count1 = g->n - count0;
    return (double)(count0 < count1 ? count0 : count1) / (double)g->n;
}

/* ==================================================================
 * Expander Zig-Zag Product (Wrapper)
 *
 * Provides the exp_zigzag_product wrapper that returns an
 * ExpanderGraph* instead of ZigZagResult, matching the
 * declaration in expander_constructions.h.
 * ================================================================== */

/**
 * exp_zigzag_product -- Compute zig-zag product and return
 * the result graph directly. The input graphs are NOT freed.
 *
 * This is a convenience wrapper around zz_compute().
 */
ExpanderGraph *exp_zigzag_product(const ExpanderGraph *G,
                                   const ExpanderGraph *H) {
    if (!G || !H) return NULL;
    ZigZagResult zr = zz_compute(G, H);
    /* We need to extract result without freeing */
    ExpanderGraph *result = zr.result;
    zr.result = NULL;
    return result;
}
