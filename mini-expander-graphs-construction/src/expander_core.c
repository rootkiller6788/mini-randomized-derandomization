/**********************************************************************
 * expander_core.c -- Core Expander Graph Operations
 *
 * Implements foundational data structures and operations for expander
 * graphs: creation, edge management, regularity checks, eigenvalue
 * computation, spectral gap analysis, Cheeger constant, expansion
 * verification, and random walk mixing time.
 *
 * Knowledge Coverage:
 *   L1: ExpanderGraph struct, GraphVector, adjacency representation
 *   L2: Spectral gap, mixing time, conductance (Cheeger constant)
 *   L3: Adjacency matrix, normalized Laplacian, eigenvalue decomposition
 *   L4: Cheeger inequality (relation between spectral and edge expansion)
 *   L5: Power iteration for eigenvalue computation, random walk simulation
 *
 * References:
 *   - Hoory, Linial, Wigderson (2006) "Expander Graphs and Their Applications"
 *   - Arora & Barak (2009) "Computational Complexity: A Modern Approach"
 *   - Alon & Spencer (2016) "The Probabilistic Method"
 **********************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "expander_core.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==================================================================
 * Graph Creation & Memory Management
 * ================================================================== */

ExpanderGraph *exp_create(size_t n, size_t d) {
    if (n == 0 || d == 0) return NULL;
    ExpanderGraph *g = (ExpanderGraph *)malloc(sizeof(ExpanderGraph));
    if (!g) return NULL;
    g->n = n;
    g->d = d;
    g->lambda = 0.0;
    g->adj = (size_t *)malloc(n * d * sizeof(size_t));
    if (!g->adj) {
        free(g);
        return NULL;
    }
    for (size_t i = 0; i < n * d; i++) {
        g->adj[i] = SIZE_MAX;  /* Sentinel: no edge */
    }
    return g;
}

void exp_free(ExpanderGraph *g) {
    if (!g) return;
    free(g->adj);
    free(g);
}

/* ==================================================================
 * Edge Management
 * ================================================================== */

bool exp_add_edge(ExpanderGraph *g, size_t u, size_t v) {
    if (!g || u >= g->n || v >= g->n || u == v) return false;

    /* Check for duplicate */
    for (size_t i = 0; i < g->d; i++) {
        if (g->adj[u * g->d + i] == v) return false;
        if (g->adj[v * g->d + i] == u) return false;
    }

    /* Find free slots */
    size_t u_slot = SIZE_MAX, v_slot = SIZE_MAX;
    for (size_t i = 0; i < g->d; i++) {
        if (g->adj[u * g->d + i] == SIZE_MAX && u_slot == SIZE_MAX) u_slot = i;
        if (g->adj[v * g->d + i] == SIZE_MAX && v_slot == SIZE_MAX) v_slot = i;
    }
    if (u_slot == SIZE_MAX || v_slot == SIZE_MAX) return false;

    g->adj[u * g->d + u_slot] = v;
    g->adj[v * g->d + v_slot] = u;
    return true;
}

/* ==================================================================
 * Regularity Check
 * ================================================================== */

bool exp_is_regular(const ExpanderGraph *g) {
    if (!g) return false;
    for (size_t u = 0; u < g->n; u++) {
        size_t deg = 0;
        for (size_t i = 0; i < g->d; i++) {
            if (g->adj[u * g->d + i] != SIZE_MAX) deg++;
        }
        if (deg != g->d) return false;
    }
    return true;
}

/* ==================================================================
 * Adjacency Matrix Construction (Internal Helpers)
 * ================================================================== */

static double **exp_adjacency_matrix(const ExpanderGraph *g) {
    if (!g) return NULL;
    size_t n = g->n, d = g->d;
    double **A = (double **)malloc(n * sizeof(double *));
    if (!A) return NULL;
    for (size_t i = 0; i < n; i++) {
        A[i] = (double *)calloc(n, sizeof(double));
        if (!A[i]) {
            for (size_t j = 0; j < i; j++) free(A[j]);
            free(A);
            return NULL;
        }
    }
    for (size_t u = 0; u < n; u++) {
        for (size_t k = 0; k < d; k++) {
            size_t v = g->adj[u * d + k];
            if (v != SIZE_MAX && v < n) A[u][v] = 1.0;
        }
    }
    return A;
}

static void exp_adjacency_matrix_free(double **A, size_t n) {
    if (!A) return;
    for (size_t i = 0; i < n; i++) free(A[i]);
    free(A);
}

/* ==================================================================
 * Power Iteration for Second Eigenvalue (Deflated)
 *
 * For d-regular graphs, the top eigenvalue is d with eigenvector
 * (1,1,...,1)/sqrt(n). We orthogonalize against this and apply
 * power iteration to extract lambda_2.
 *
 * Reference: Golub & Van Loan (2013) "Matrix Computations", Sec 8.2
 * ================================================================== */

static double exp_power_iteration_deflated(const ExpanderGraph *g, int iters) {
    if (!g || g->n <= 1) return 0.0;
    size_t n = g->n;
    double **A = exp_adjacency_matrix(g);
    if (!A) return 0.0;

    double *v = (double *)malloc(n * sizeof(double));
    double *w = (double *)malloc(n * sizeof(double));
    if (!v || !w) {
        free(v); free(w);
        exp_adjacency_matrix_free(A, n);
        return 0.0;
    }

    /* Initialize random vector orthogonal to all-ones */
    for (size_t i = 0; i < n; i++)
        v[i] = ((double)rand() / RAND_MAX) - 0.5;

    /* Orthogonalize against all-ones */
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) sum += v[i];
    sum /= (double)n;
    for (size_t i = 0; i < n; i++) v[i] -= sum;

    /* Normalize */
    double norm = 0.0;
    for (size_t i = 0; i < n; i++) norm += v[i] * v[i];
    norm = sqrt(norm);
    if (norm > 1e-15)
        for (size_t i = 0; i < n; i++) v[i] /= norm;

    double lambda = 0.0;
    for (int t = 0; t < iters; t++) {
        /* w = A * v */
        for (size_t i = 0; i < n; i++) {
            w[i] = 0.0;
            for (size_t j = 0; j < n; j++)
                w[i] += A[i][j] * v[j];
        }
        /* Deflate: remove all-ones component */
        double sum_w = 0.0;
        for (size_t i = 0; i < n; i++) sum_w += w[i];
        double proj = sum_w / (double)n;
        for (size_t i = 0; i < n; i++) w[i] -= proj;

        /* Normalize */
        norm = 0.0;
        for (size_t i = 0; i < n; i++) norm += w[i] * w[i];
        norm = sqrt(norm);
        if (norm < 1e-15) break;
        for (size_t i = 0; i < n; i++) v[i] = w[i] / norm;

        /* Rayleigh quotient */
        double num = 0.0;
        for (size_t i = 0; i < n; i++) {
            double Av = 0.0;
            for (size_t j = 0; j < n; j++) Av += A[i][j] * v[j];
            num += v[i] * Av;
        }
        lambda = num;  /* v is normalized */
    }

    free(v); free(w);
    exp_adjacency_matrix_free(A, n);
    return lambda;
}

/* ==================================================================
 * Eigenvalue and Spectral Analysis
 * ================================================================== */

double exp_second_eigenvalue(const ExpanderGraph *g) {
    if (!g || g->n <= 1) return 0.0;
    return exp_power_iteration_deflated(g, 200);
}

double exp_spectral_gap(const ExpanderGraph *g) {
    if (!g) return 0.0;
    double lam2 = exp_second_eigenvalue(g);
    double gap = (double)g->d - lam2;
    return gap > 0 ? gap : 0.0;
}

/* ==================================================================
 * Expander Verification
 * ================================================================== */

bool exp_is_expander(const ExpanderGraph *g, double gamma) {
    if (!g || gamma <= 0.0 || gamma > 1.0) return false;
    if (g->n <= 1) return false;
    double lam2 = exp_second_eigenvalue(g);
    return lam2 <= (double)g->d * (1.0 - gamma);
}

bool exp_is_lambda_expander(const ExpanderGraph *g, double lambda) {
    if (!g || lambda <= 0.0) return false;
    double lam2 = exp_second_eigenvalue(g);
    return lam2 <= lambda;
}

/* ==================================================================
 * Cheeger Constant (Edge Expansion / Conductance)
 *
 * h(G) = min_{S,|S|<=n/2} |E(S,V\S)| / |S|
 *
 * Spectral bounds via Cheeger inequality (Alon-Milman 1985):
 *   gamma/2 <= h(G) <= sqrt(2*d*gamma)
 * where gamma = d - lambda_2 is the spectral gap.
 * ================================================================== */

double exp_cheeger_constant(const ExpanderGraph *g) {
    if (!g || g->n <= 1) return 0.0;
    double lam2 = exp_second_eigenvalue(g);
    double gap = (double)g->d - lam2;
    return sqrt(2.0 * (double)g->d * gap);
}

double exp_cheeger_lower_bound(const ExpanderGraph *g) {
    if (!g) return 0.0;
    double lam2 = exp_second_eigenvalue(g);
    double gap = (double)g->d - lam2;
    return gap / 2.0;
}

void exp_cheeger_bounds(const ExpanderGraph *g, double *lower, double *upper) {
    if (!g || !lower || !upper) return;
    double lam2 = exp_second_eigenvalue(g);
    double gap = (double)g->d - lam2;
    *lower = (gap > 0) ? gap / 2.0 : 0.0;
    *upper = (gap > 0) ? sqrt(2.0 * (double)g->d * gap) : 0.0;
}

/* ==================================================================
 * Expansion Verification (Vertex Expansion)
 * ================================================================== */

bool exp_verify_expansion(const ExpanderGraph *g, size_t S, double phi) {
    if (!g || S >= g->n || phi < 0) return false;
    size_t max_n = g->n;
    char *in_S = (char *)calloc(max_n, sizeof(char));
    char *in_N = (char *)calloc(max_n, sizeof(char));
    if (!in_S || !in_N) { free(in_S); free(in_N); return false; }

    in_S[S] = 1;
    for (size_t k = 0; k < g->d; k++) {
        size_t v = g->adj[S * g->d + k];
        if (v != SIZE_MAX && v < max_n) in_N[v] = 1;
    }

    size_t exclusive_n = 0;
    for (size_t v = 0; v < max_n; v++)
        if (in_N[v] && !in_S[v]) exclusive_n++;

    free(in_S); free(in_N);
    return (double)exclusive_n >= phi;
}

double exp_vertex_expansion_ratio(const ExpanderGraph *g, const char *subset) {
    if (!g || !subset) return -1.0;
    char *Nset = (char *)calloc(g->n, sizeof(char));
    if (!Nset) return -1.0;

    size_t s_size = 0;
    for (size_t i = 0; i < g->n; i++) {
        if (subset[i]) {
            Nset[i] = 1;
            s_size++;
            for (size_t k = 0; k < g->d; k++) {
                size_t v = g->adj[i * g->d + k];
                if (v != SIZE_MAX && v < g->n) Nset[v] = 1;
            }
        }
    }
    if (s_size == 0) { free(Nset); return -1.0; }

    size_t n_size = 0;
    for (size_t i = 0; i < g->n; i++) if (Nset[i]) n_size++;
    free(Nset);
    return (double)(n_size - s_size) / (double)s_size;
}

/* ==================================================================
 * Random Walk Mixing Time
 *
 * tau(eps) <= (1/(1-lambda)) * log(n/eps)
 * where lambda = lambda_2 / d is the normalized second eigenvalue.
 *
 * Reference: Levin, Peres, Wilmer (2017) "Markov Chains and Mixing Times"
 * ================================================================== */

size_t exp_random_walk_mix_time(const ExpanderGraph *g, double eps) {
    if (!g || eps <= 0.0 || eps >= 1.0) return 0;
    if (g->n <= 1) return 0;
    double lam2 = exp_second_eigenvalue(g);
    double lam = lam2 / (double)g->d;
    if (lam >= 1.0) lam = 0.999;
    double gap = 1.0 - lam;
    double tau = log((double)g->n / eps) / gap;
    return (size_t)ceil(tau);
}

size_t exp_mixing_time_precise(const ExpanderGraph *g, double eps) {
    if (!g || eps <= 0.0) return 0;
    if (g->n <= 1) return 0;
    double lam2 = exp_second_eigenvalue(g);
    double lam = lam2 / (double)g->d;
    if (lam >= 1.0) lam = 0.999;
    double gap = 1.0 - lam;
    double min_pi = 1.0 / (double)g->n;
    double tau = log(1.0 / (eps * sqrt(min_pi))) / gap;
    return (size_t)ceil(tau);
}

/* ==================================================================
 * Random Walk Simulation
 * ================================================================== */

size_t exp_random_walk_step(const ExpanderGraph *g, size_t current) {
    if (!g || current >= g->n) return current;
    if (rand() < RAND_MAX / 2) return current;
    int idx = rand() % (int)g->d;
    size_t next = g->adj[current * g->d + idx];
    if (next == SIZE_MAX || next >= g->n) return current;
    return next;
}

void exp_random_walk_distribution(const ExpanderGraph *g, size_t start,
                                   size_t steps, size_t samples, double *dist) {
    if (!g || !dist || start >= g->n) return;
    memset(dist, 0, g->n * sizeof(double));
    for (size_t s = 0; s < samples; s++) {
        size_t pos = start;
        for (size_t t = 0; t < steps; t++)
            pos = exp_random_walk_step(g, pos);
        dist[pos] += 1.0;
    }
    for (size_t i = 0; i < g->n; i++) dist[i] /= (double)samples;
}

double exp_total_variation_distance(const double *p, size_t n) {
    if (!p || n == 0) return 1.0;
    double tv = 0.0;
    for (size_t i = 0; i < n; i++)
        tv += fabs(p[i] - 1.0 / (double)n);
    return 0.5 * tv;
}

/* ==================================================================
 * Graph Statistics
 * ================================================================== */

size_t exp_num_edges(const ExpanderGraph *g) {
    if (!g) return 0;
    size_t edge_count = 0;
    for (size_t u = 0; u < g->n; u++)
        for (size_t k = 0; k < g->d; k++)
            if (g->adj[u * g->d + k] != SIZE_MAX) edge_count++;
    return edge_count / 2;
}

size_t exp_diameter(const ExpanderGraph *g) {
    if (!g || g->n <= 1) return 0;
    double lam2 = exp_second_eigenvalue(g);
    if (lam2 >= (double)g->d || lam2 <= 0) return g->n;
    double ratio = (double)g->d / lam2;
    return (size_t)ceil(log((double)(g->n - 1)) / log(ratio));
}

size_t exp_girth_lower_bound(const ExpanderGraph *g) {
    if (!g || g->n < 3) return 0;
    double lam2 = exp_second_eigenvalue(g);
    double d = (double)g->d;
    if (lam2 >= 2.0 * sqrt(d - 1.0)) return 3;
    double ratio = lam2 / (d - 1.0);
    if (ratio <= 0) return 3;
    return (size_t)floor(log((double)g->n) / log(1.0 / ratio));
}

/* ==================================================================
 * Exact Edge Expansion (Small Graphs Only)
 * ================================================================== */

double exp_edge_expansion_exact(const ExpanderGraph *g) {
    if (!g || g->n > 20 || g->n == 0) return -1.0;
    double best = INFINITY;
    size_t n = g->n, half = n / 2;

    for (size_t mask = 1; mask < (1ULL << n); mask++) {
        size_t s_size = 0;
        for (size_t i = 0; i < n; i++)
            if (mask & (1ULL << i)) s_size++;
        if (s_size > half || s_size == 0) continue;

        size_t crossing = 0;
        for (size_t u = 0; u < n; u++) {
            if (!(mask & (1ULL << u))) continue;
            for (size_t k = 0; k < g->d; k++) {
                size_t v = g->adj[u * g->d + k];
                if (v != SIZE_MAX && v < n && !(mask & (1ULL << v)))
                    crossing++;
            }
        }
        double ratio = (double)crossing / (double)s_size;
        if (ratio < best) best = ratio;
    }
    return (best == INFINITY) ? 0.0 : best;
}

/* ==================================================================
 * GraphVector Operations
 * ================================================================== */

GraphVector *gv_create(size_t n) {
    if (n == 0) return NULL;
    GraphVector *v = (GraphVector *)malloc(sizeof(GraphVector));
    if (!v) return NULL;
    v->n = n;
    v->vec = (double *)calloc(n, sizeof(double));
    if (!v->vec) { free(v); return NULL; }
    return v;
}

void gv_free(GraphVector *v) {
    if (!v) return;
    free(v->vec);
    free(v);
}

double gv_norm(const GraphVector *v) {
    if (!v) return 0.0;
    double sum = 0.0;
    for (size_t i = 0; i < v->n; i++) sum += v->vec[i] * v->vec[i];
    return sqrt(sum);
}

void gv_normalize(GraphVector *v) {
    if (!v) return;
    double nrm = gv_norm(v);
    if (nrm < 1e-15) return;
    for (size_t i = 0; i < v->n; i++) v->vec[i] /= nrm;
}

double gv_dot(const GraphVector *a, const GraphVector *b) {
    if (!a || !b || a->n != b->n) return 0.0;
    double sum = 0.0;
    for (size_t i = 0; i < a->n; i++) sum += a->vec[i] * b->vec[i];
    return sum;
}

void exp_apply_adjacency(const ExpanderGraph *g, const GraphVector *v,
                          GraphVector *result) {
    if (!g || !v || !result || v->n != g->n || result->n != g->n) return;
    for (size_t i = 0; i < g->n; i++) {
        result->vec[i] = 0.0;
        for (size_t k = 0; k < g->d; k++) {
            size_t j = g->adj[i * g->d + k];
            if (j != SIZE_MAX && j < g->n)
                result->vec[i] += v->vec[j];
        }
    }
}

double exp_rayleigh_quotient(const ExpanderGraph *g, const GraphVector *v) {
    if (!g || !v || v->n != g->n) return 0.0;
    GraphVector *Av = gv_create(g->n);
    if (!Av) return 0.0;
    exp_apply_adjacency(g, v, Av);
    double num = gv_dot(v, Av);
    double den = gv_dot(v, v);
    gv_free(Av);
    if (den < 1e-15) return 0.0;
    return num / den;
}

/* ==================================================================
 * Advanced Spectral Quantities
 * ================================================================== */

double exp_collision_probability(const ExpanderGraph *g, size_t steps) {
    if (!g || g->n == 0) return 1.0;
    double *dist = (double *)calloc(g->n, sizeof(double));
    if (!dist) return 1.0;
    dist[0] = 1.0;

    for (size_t t = 0; t < steps; t++) {
        double *next = (double *)calloc(g->n, sizeof(double));
        if (!next) { free(dist); return 1.0; }
        for (size_t u = 0; u < g->n; u++) {
            double stay = 0.5 * dist[u];
            double move = 0.5 * dist[u] / (double)g->d;
            next[u] += stay;
            for (size_t k = 0; k < g->d; k++) {
                size_t v = g->adj[u * g->d + k];
                if (v != SIZE_MAX && v < g->n) next[v] += move;
            }
        }
        free(dist);
        dist = next;
    }

    double cp = 0.0;
    for (size_t i = 0; i < g->n; i++) cp += dist[i] * dist[i];
    free(dist);
    return cp;
}

double exp_entropy_rate(const ExpanderGraph *g) {
    if (!g || g->d == 0) return 0.0;
    double d = (double)g->d;
    return 0.5 + 0.5 * log2(2.0 * d);
}

/* ==================================================================
 * Graph Property Queries
 * ================================================================== */

size_t exp_degree(const ExpanderGraph *g, size_t u) {
    if (!g || u >= g->n) return 0;
    size_t deg = 0;
    for (size_t i = 0; i < g->d; i++)
        if (g->adj[u * g->d + i] != SIZE_MAX) deg++;
    return deg;
}

bool exp_has_edge(const ExpanderGraph *g, size_t u, size_t v) {
    if (!g || u >= g->n || v >= g->n) return false;
    for (size_t i = 0; i < g->d; i++)
        if (g->adj[u * g->d + i] == v) return true;
    return false;
}

size_t exp_neighbors(const ExpanderGraph *g, size_t u, size_t *neighbors) {
    if (!g || !neighbors || u >= g->n) return 0;
    size_t count = 0;
    for (size_t i = 0; i < g->d; i++) {
        size_t v = g->adj[u * g->d + i];
        if (v != SIZE_MAX && v < g->n) neighbors[count++] = v;
    }
    return count;
}

bool exp_is_bipartite(const ExpanderGraph *g) {
    if (!g || g->n == 0) return false;
    char *color = (char *)malloc(g->n);
    if (!color) return false;
    for (size_t i = 0; i < g->n; i++) color[i] = (char)-1;

    for (size_t start = 0; start < g->n; start++) {
        if (color[start] != (char)-1) continue;
        color[start] = 0;
        size_t *queue = (size_t *)malloc(g->n * sizeof(size_t));
        if (!queue) { free(color); return false; }
        size_t head = 0, tail = 0;
        queue[tail++] = start;

        while (head < tail) {
            size_t u = queue[head++];
            char next_color = (char)(1 - color[u]);
            for (size_t k = 0; k < g->d; k++) {
                size_t v = g->adj[u * g->d + k];
                if (v == SIZE_MAX || v >= g->n) continue;
                if (color[v] == (char)-1) {
                    color[v] = next_color;
                    queue[tail++] = v;
                } else if (color[v] != next_color) {
                    free(queue); free(color);
                    return false;
                }
            }
        }
        free(queue);
    }
    free(color);
    return true;
}

bool exp_is_connected(const ExpanderGraph *g) {
    if (!g || g->n == 0) return false;
    char *visited = (char *)calloc(g->n, sizeof(char));
    if (!visited) return false;

    size_t *queue = (size_t *)malloc(g->n * sizeof(size_t));
    if (!queue) { free(visited); return false; }
    size_t head = 0, tail = 0;
    visited[0] = 1;
    queue[tail++] = 0;

    while (head < tail) {
        size_t u = queue[head++];
        for (size_t k = 0; k < g->d; k++) {
            size_t v = g->adj[u * g->d + k];
            if (v != SIZE_MAX && v < g->n && !visited[v]) {
                visited[v] = 1;
                queue[tail++] = v;
            }
        }
    }

    bool connected = (tail == g->n);
    free(queue); free(visited);
    return connected;
}

void exp_print_summary(const ExpanderGraph *g, FILE *fp) {
    if (!g || !fp) return;
    fprintf(fp, "ExpanderGraph: n=%zu, d=%zu\n", g->n, g->d);
    fprintf(fp, "  Edges: %zu\n", exp_num_edges(g));
    if (exp_is_regular(g)) {
        double lam2 = exp_second_eigenvalue(g);
        double gap = exp_spectral_gap(g);
        double h_lo, h_hi;
        exp_cheeger_bounds(g, &h_lo, &h_hi);
        fprintf(fp, "  Regular: yes\n");
        fprintf(fp, "  lambda_2: %.6f\n", lam2);
        fprintf(fp, "  Spectral gap: %.6f\n", gap);
        fprintf(fp, "  Cheeger bounds: [%.6f, %.6f]\n", h_lo, h_hi);
        fprintf(fp, "  Diameter bound: %zu\n", exp_diameter(g));
        fprintf(fp, "  Mixing time (eps=0.01): %zu\n",
                exp_random_walk_mix_time(g, 0.01));
        fprintf(fp, "  Connected: %s\n",
                exp_is_connected(g) ? "yes" : "no");
        fprintf(fp, "  Bipartite: %s\n",
                exp_is_bipartite(g) ? "yes" : "no");
        double alon = exp_alon_boppana_bound_static(g->d);
        fprintf(fp, "  Alon-Boppana bound: %.6f\n", alon);
        if (lam2 <= 2.0 * sqrt((double)g->d - 1.0))
            fprintf(fp, "  Ramanujan: yes\n");
        else
            fprintf(fp, "  Ramanujan: no\n");
    } else {
        fprintf(fp, "  Regular: no\n");
    }
}

double exp_alon_boppana_bound_static(size_t d) {
    if (d <= 1) return (double)d;
    return 2.0 * sqrt((double)(d - 1));
}

double exp_compare_to_optimal(const ExpanderGraph *g) {
    if (!g || g->d <= 1) return INFINITY;
    double lam2 = exp_second_eigenvalue(g);
    double optimal = 2.0 * sqrt((double)(g->d - 1));
    if (optimal == 0.0) return INFINITY;
    return lam2 / optimal;
}
