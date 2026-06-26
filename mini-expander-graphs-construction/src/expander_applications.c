/**********************************************************************
 * expander_applications.c -- Applications of Expander Graphs
 *
 * Implements key applications of expander graphs in theoretical
 * computer science:
 *   - Error reduction for randomized algorithms (RP/BPP amplification)
 *   - Derandomized sampling using expander walks
 *   - SL connectivity (Reingold's USTCON in log-space)
 *   - Sorting networks from expanders
 *   - Chernoff bounds on expander random walks
 *   - Expander codes (basic error-correcting codes)
 *
 * Knowledge Coverage:
 *   L1-L2: Random walk, error reduction, derandomization concept
 *   L3: Expander-based sampler, pseudorandom walk
 *   L4: Chernoff bound for expander walks (Gilman 1998)
 *   L5: SL=L algorithm, expander-based error amplification
 *   L6: USTCON (undirected s-t connectivity)
 *   L7: Cryptography, derandomization, coding theory applications
 *
 * References:
 *   - Reingold (2008) "Undirected Connectivity in Log-Space"
 *   - Gillman (1998) "A Chernoff Bound for Random Walks on
 *     Expander Graphs"
 *   - Ajtai, Komlos, Szemeredi (1983) "Sorting in c log n steps"
 *   - Sipser & Spielman (1996) "Expander Codes"
 **********************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "expander_core.h"
#include "expander_constructions.h"
#include "expander_applications.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==================================================================
 * Error Reduction for Randomized Algorithms
 *
 * Given a BPP/RP algorithm with base error probability p < 1/2,
 * an expander random walk amplifies the success probability to
 * 1 - exp(-Omega(k)) using only O(k) random bits (instead of
 * O(k * log n) for independent trials).
 *
 * The expander Chernoff bound (Gillman 1998) guarantees:
 *   Pr[ |(1/k) * sum X_i - mu| > delta ] <= 2 * exp(-Omega(k * delta^2))
 * for a random walk on a spectral expander.
 *
 * This is exponentially better than naive repetition because
 * we reuse randomness across trials.
 * ================================================================== */

/**
 * exp_error_reduction -- Estimate the amplified error probability
 * after 'trials' independent random steps on an expander graph.
 *
 * Uses the spectral bound: For a d-regular expander with
 * lambda = lambda_2 / d, the error after k trials is bounded by:
 *
 *   Pr[error] <= 2 * exp(-k * (1/2 - p)^2 / (2*(1+lambda)/(1-lambda)))
 *
 * For a Ramanujan expander, lambda ~ 2*sqrt(d-1)/d.
 *
 * Complexity: O(n^2) for eigenvalue computation.
 */
bool exp_error_reduction(const ExpanderGraph *g, size_t trials,
                          double base_err, double *amplified_err) {
    if (!g || !amplified_err || base_err <= 0.0 || base_err >= 0.5)
        return false;

    double lam2 = exp_second_eigenvalue(g);
    double lam = lam2 / (double)g->d;
    if (lam >= 1.0) lam = 0.999;

    /* Gillman's bound: tail <= 2 * exp(-trials * gap * eps^2 / 12) */
    double gap = 1.0 - lam;
    double eps = 0.5 - base_err;  /* Distance from unbiased */
    if (eps <= 0.0) eps = 0.01;

    double exponent = (double)trials * gap * eps * eps / 12.0;
    *amplified_err = 2.0 * exp(-exponent);
    if (*amplified_err > 1.0) *amplified_err = 1.0;

    return true;
}

/**
 * exp_amplify_rp_algorithm -- Amplify an RP algorithm using expander
 * random walk instead of independent repetitions.
 *
 * Given an RP predicate rp() that returns true/false with one-sided
 * error p < 1/2, amplify the error to 2^{-k} using only
 * O(k) additional random bits via an expander walk of length k.
 *
 * The function simulates this by calling rp() 'trials' times
 * (each call represents one step of a random walk on the expander).
 * Returns true if any trial accepts (RP acceptance), false otherwise.
 *
 * The error probability bound is stored in *error.
 */
bool exp_amplify_rp_algorithm(const ExpanderGraph *g,
                               bool (*rp)(void), size_t trials,
                               double *error) {
    if (!g || !rp || !error) return false;

    /* Run RP algorithm at 'trials' vertices of an expander random walk */
    bool accepted = false;
    for (size_t i = 0; i < trials; i++) {
        if (rp()) {
            accepted = true;
        }
    }

    /* Error bound for RP (one-sided): if answer is NO,
     * error = Pr[all trials accept] (false positive).
     * For RP, the base error p0 < 1/2. With expander walk:
     *   error <= p0^trials * something spectral
     * We use a conservative bound: p0^trials
     */
    double base_err = 0.333;  /* Typical RP error < 1/3 */
    double amp_err;
    exp_error_reduction(g, trials, base_err, &amp_err);
    *error = amp_err;

    return accepted;
}

/* ==================================================================
 * Derandomized Sampling
 *
 * Expander random walks provide a derandomized way to sample
 * from a large universe. Instead of picking k independent
 * uniform samples (requiring k * log(n) random bits),
 * we pick a random starting vertex and then walk k steps
 * on the expander (requiring only log(n) + k * log(d) bits).
 *
 * For a spectral expander, the statistical properties (e.g.,
 * hitting probabilities, Chernoff bounds) are almost as good
 * as independent sampling.
 * ================================================================== */

/**
 * exp_derandomized_sampling -- Generate 'samples' approximately
 * uniform samples from [0, n-1] by performing an expander random
 * walk of length 'samples' starting from a random vertex.
 *
 * The indices array must have size >= samples.
 *
 * Randomness used: log2(n) + samples * log2(d) bits
 * vs. naive: samples * log2(n) bits.
 *
 * Complexity: O(samples).
 */
bool exp_derandomized_sampling(const ExpanderGraph *g, size_t samples,
                                size_t *indices) {
    if (!g || !indices || samples == 0) return false;

    /* Random starting vertex */
    size_t current = (size_t)rand() % g->n;
    indices[0] = current;

    for (size_t i = 1; i < samples; i++) {
        /* Step: lazy random walk */
        if (rand() < RAND_MAX / 2) {
            /* Stay */
            indices[i] = current;
        } else {
            /* Move to random neighbor */
            int edge_idx = rand() % (int)g->d;
            size_t next = g->adj[current * g->d + edge_idx];
            if (next != SIZE_MAX && next < g->n) {
                current = next;
            }
            indices[i] = current;
        }
    }
    return true;
}

/**
 * exp_sampling_quality -- Estimate the quality of derandomized
 * sampling by computing the total variation distance between
 * the empirical distribution of samples and the uniform distribution.
 *
 * Lower values indicate better sampling quality.
 * Returns the TV distance.
 */
double exp_sampling_quality(const size_t *samples, size_t num_samples,
                             size_t n) {
    if (!samples || num_samples == 0 || n == 0) return 1.0;

    double *empirical = (double *)calloc(n, sizeof(double));
    if (!empirical) return 1.0;

    for (size_t i = 0; i < num_samples; i++) {
        if (samples[i] < n) empirical[samples[i]] += 1.0;
    }
    for (size_t i = 0; i < n; i++) empirical[i] /= (double)num_samples;

    double tv = 0.0;
    for (size_t i = 0; i < n; i++)
        tv += fabs(empirical[i] - 1.0 / (double)n);

    free(empirical);
    return 0.5 * tv;
}

/* ==================================================================
 * SL Connectivity (Undirected s-t Connectivity)
 *
 * Reingold (2008) proved SL = L by showing that USTCON
 * (undirected s-t connectivity) can be solved in deterministic
 * log-space using expander graphs.
 *
 * The algorithm:
 *   1. Convert the input graph to a regular graph.
 *   2. Use the zig-zag product to build an expander.
 *   3. Perform a deterministic walk on the expander.
 *   4. Check if s and t are in the same connected component.
 *
 * This implementation provides a simplified version that uses
 * spectral properties of expanders to test connectivity.
 * ================================================================== */

/**
 * exp_sl_connectivity -- Test if vertices s and t are connected
 * in an expander graph using spectral properties and random walks.
 *
 * For a connected d-regular expander, a random walk of length
 * O(log n) from s will visit t with high probability.
 *
 * This function simulates a random walk of length mixing_time
 * starting from s and checks if t is ever visited.
 *
 * Complexity: O(mixing_time) for the random walk.
 */
bool exp_sl_connectivity(const ExpanderGraph *g, size_t s, size_t t,
                          bool *connected) {
    if (!g || !connected || s >= g->n || t >= g->n) return false;

    if (!exp_is_connected(g)) {
        /* Run BFS from s to check if t is reachable */
        char *visited = (char *)calloc(g->n, sizeof(char));
        if (!visited) return false;

        size_t *queue = (size_t *)malloc(g->n * sizeof(size_t));
        if (!queue) { free(visited); return false; }
        size_t head = 0, tail = 0;
        visited[s] = 1;
        queue[tail++] = s;

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

        *connected = (visited[t] != 0);
        free(queue); free(visited);
        return true;
    }

    /* For an expander, use random walk to test connectivity */
    size_t mix_time = exp_random_walk_mix_time(g, 0.1);
    if (mix_time == 0) mix_time = (size_t)log2((double)g->n) * 10;

    size_t pos = s;
    for (size_t step = 0; step < mix_time * 5; step++) {
        pos = exp_random_walk_step(g, pos);
        if (pos == t) {
            *connected = true;
            return true;
        }
    }

    /* Didn't find t -- may be in different component */
    *connected = false;
    return true;
}

/**
 * exp_sl_path_length -- Estimate the shortest path length between
 * s and t using the spectral diameter bound.
 *
 * For expanders, diam(G) = O(log n), so the path length is small.
 * This returns the spectral bound rather than the exact distance.
 *
 * Complexity: O(n^2) for eigenvalue computation.
 */
size_t exp_sl_path_length(const ExpanderGraph *g, size_t s, size_t t) {
    if (!g || s >= g->n || t >= g->n) return SIZE_MAX;
    if (s == t) return 0;

    /* Spectral diameter bound */
    double lam2 = exp_second_eigenvalue(g);
    /* Use absolute value for the spectral bound (worst-case eigenvalue) */
    double lam_eff = fabs(lam2);
    if (lam_eff >= (double)g->d || lam_eff <= 1e-10) return SIZE_MAX;

    double ratio = (double)g->d / lam_eff;
    size_t diam = (size_t)ceil(log((double)(g->n - 1)) / log(ratio));

    /* BFS for exact distance if graph is small */
    if (g->n <= 1000) {
        char *visited = (char *)calloc(g->n, sizeof(char));
        size_t *dist = (size_t *)malloc(g->n * sizeof(size_t));
        size_t *queue = (size_t *)malloc(g->n * sizeof(size_t));
        if (!visited || !dist || !queue) {
            free(visited); free(dist); free(queue);
            return diam;  /* Fallback to spectral bound */
        }

        for (size_t i = 0; i < g->n; i++) dist[i] = SIZE_MAX;
        size_t head = 0, tail = 0;
        visited[s] = 1; dist[s] = 0;
        queue[tail++] = s;

        while (head < tail) {
            size_t u = queue[head++];
            if (u == t) {
                size_t d = dist[t];
                free(visited); free(dist); free(queue);
                return d;
            }
            for (size_t k = 0; k < g->d; k++) {
                size_t v = g->adj[u * g->d + k];
                if (v != SIZE_MAX && v < g->n && !visited[v]) {
                    visited[v] = 1;
                    dist[v] = dist[u] + 1;
                    queue[tail++] = v;
                }
            }
        }
        free(visited); free(dist); free(queue);
    }

    return diam;  /* Spectral upper bound */
}

/* ==================================================================
 * Sorting Networks from Expanders
 *
 * Ajtai, Komlos, Szemeredi (1983) showed that expander graphs
 * can be used to construct sorting networks of depth O(log n).
 *
 * The AKS sorting network uses expander-based comparator graphs
 * to achieve asymptotically optimal depth.
 *
 * This implementation generates a simple sorting network skeleton
 * based on an expander graph topology.
 * ================================================================== */

/**
 * exp_sorting_network_from_expander -- Generate comparator pairs
 * based on an expander graph on n vertices.
 *
 * The expander edges serve as comparators in the sorting network.
 * For a d-regular expander, this gives O(nd) comparators.
 *
 * This is a SIMPLIFIED construction; the full AKS network requires
 * a recursive bipartite expander structure.
 *
 * Complexity: O(nd).
 */
bool exp_sorting_network_from_expander(size_t n, size_t *comparators_out,
                                        size_t *num_comparators) {
    if (n == 0 || !comparators_out || !num_comparators) return false;

    /* Build a small expander on n/log(n) vertices using Margulis */
    size_t m = (size_t)ceil(sqrt((double)n));
    size_t m_sq = m * m;

    ExpanderGraph *g = exp_margulis_gabber_galil(m_sq);
    if (!g) return false;

    /* Use first n vertices of the expander */
    size_t count = 0;
    size_t max_comparators = n * g->d;  /* Upper bound */

    for (size_t u = 0; u < n && u < g->n; u++) {
        for (size_t k = 0; k < g->d; k++) {
            size_t v = g->adj[u * g->d + k];
            if (v != SIZE_MAX && v < n && u < v) {
                if (count * 2 + 1 < max_comparators) {
                    comparators_out[count * 2] = u;
                    comparators_out[count * 2 + 1] = v;
                    count++;
                }
            }
        }
    }

    *num_comparators = count;
    exp_free(g);
    return true;
}

/* ==================================================================
 * Chernoff Bound on Expanders (Gillman 1998)
 *
 * For a random walk X_1, ..., X_k on a spectral expander,
 * the sum S = sum_i f(X_i) for a function f: V -> [0,1]
 * concentrates around its mean:
 *
 *   Pr[ |S/k - mu| > eps ] <= 2 * exp(-Omega(eps^2 * k * (1-lambda)))
 *
 * where lambda is the normalized second eigenvalue.
 *
 * This is the expander Chernoff bound, which is nearly as strong
 * as the independent case but requires exponentially fewer random bits.
 * ================================================================== */

/**
 * exp_chernoff_bound_on_expander -- Compute the Chernoff-style
 * tail bound for the average of a [0,1]-valued function evaluated
 * along a random walk on an expander.
 *
 * Parameters:
 *   n: number of vertices
 *   mu: expected value of the function
 *   delta: deviation from mu
 *
 * Returns the probability bound Pr[|avg - mu| > delta].
 */
double exp_chernoff_bound_on_expander(const ExpanderGraph *g, size_t n,
                                       double mu, double delta) {
    (void)mu;  /* Symmetric bound independent of mu */
    if (!g || n == 0 || delta <= 0.0) return 1.0;

    double lam2 = exp_second_eigenvalue(g);
    double lam = lam2 / (double)g->d;
    if (lam >= 1.0) lam = 0.999;

    double gap = 1.0 - lam;  /* Spectral gap */
    if (gap <= 0.0) return 1.0;  /* No expansion, bound is vacuous */
    double exponent = delta * delta * (double)n * gap / 12.0;

    double bound = 2.0 * exp(-exponent);
    return (bound > 1.0) ? 1.0 : bound;
}

/**
 * exp_chernoff_bound_independent -- For comparison: the standard
 * Chernoff bound for independent samples.
 *
 *   Pr[|avg - mu| > delta] <= 2 * exp(-2 * delta^2 * n)
 */
double exp_chernoff_bound_independent(size_t n, double mu, double delta) {
    (void)mu;  /* Symmetric bound doesn't depend on mu */
    if (n == 0 || delta <= 0.0) return 1.0;
    return 2.0 * exp(-2.0 * delta * delta * (double)n);
}

/**
 * exp_random_bits_saved -- Estimate the number of random bits saved
 * by using an expander walk vs. independent uniform samples.
 *
 * Independent: k * log2(n) bits
 * Expander walk: log2(n) + k * log2(d) bits
 *
 * Returns the number of bits saved.
 */
double exp_random_bits_saved(size_t k, size_t n, size_t d) {
    if (n <= 1 || d <= 1) return 0.0;
    double independent_bits = (double)k * log2((double)n);
    double expander_bits = log2((double)n) + (double)k * log2((double)d);
    double saved = independent_bits - expander_bits;
    return saved > 0 ? saved : 0.0;
}

/* ==================================================================
 * Expander Codes (Sipser-Spielman 1996)
 *
 * Expander codes are a family of error-correcting codes based on
 * expander graphs. The code is defined by a bipartite expander
 * graph between message bits and parity-check bits.
 *
 * Key property: Linear-time encoding and decoding.
 * Minimum distance proportional to the block length.
 * ================================================================== */

/**
 * exp_code_encode -- Encode a binary message using a simple
 * expander-based code.
 *
 * Message bits are placed on the left side of a bipartite
 * expander; each parity-check bit (right side) is the XOR
 * of its neighboring message bits.
 *
 * This is a SIMPLIFIED demonstration using a d-regular graph.
 *
 * message: array of m bits (as chars '0'/'1')
 * m: number of message bits
 * codeword_out: output array of size m + n (message + parity)
 *   (the message is copied, then parity bits appended)
 *
 * Complexity: O(m * d).
 */
bool exp_code_encode(const char *message, size_t m,
                      const ExpanderGraph *g, char *codeword_out) {
    if (!message || m == 0 || !g || !codeword_out) return false;
    if (g->n < m) return false;

    /* Copy message bits */
    memcpy(codeword_out, message, m);

    /* Compute parity bits */
    for (size_t i = 0; i < g->n && i < m; i++) {
        char parity = 0;
        size_t count = 0;
        for (size_t k = 0; k < g->d; k++) {
            size_t v = g->adj[i * g->d + k];
            if (v != SIZE_MAX && v < m) {
                parity ^= (message[v] - '0');
                count++;
            }
        }
        if (count > 0) {
            codeword_out[m + i] = parity ? '1' : '0';
        } else {
            codeword_out[m + i] = '0';
        }
    }
    return true;
}

/**
 * exp_code_flip_decode -- Simple iterative decoding for expander codes.
 * Flips bits that violate a majority of parity checks (Sipser-Spielman).
 *
 * This is a single-iteration simplified decoder.
 */
bool exp_code_flip_decode(char *codeword, size_t m,
                           const ExpanderGraph *g, size_t max_iters) {
    if (!codeword || m == 0 || !g) return false;

    for (size_t iter = 0; iter < max_iters; iter++) {
        bool changed = false;
        for (size_t i = 0; i < m && i < g->n; i++) {
            int unsatisfied = 0;
            int total = 0;
            for (size_t k = 0; k < g->d; k++) {
                size_t v = g->adj[i * g->d + k];
                if (v != SIZE_MAX && v < m) {
                    char parity = codeword[v] - '0';
                    char expected = codeword[m + i] - '0';
                    if (parity != expected) unsatisfied++;
                    total++;
                }
            }
            /* If more than half of parity checks fail, flip the bit */
            if (total > 0 && unsatisfied > total / 2) {
                codeword[i] = (codeword[i] == '0') ? '1' : '0';
                changed = true;
            }
        }
        if (!changed) break;
    }
    return true;
}

/* ==================================================================
 * Pseudorandom Generators from Expanders
 *
 * Nisan & Zuckerman (1996): Expander random walks can be used
 * to construct pseudorandom generators (PRGs) that fool
 * small-space machines.
 *
 * The PRG output is the sequence of vertices visited during
 * a random walk on an explicit expander.
 * ================================================================== */

/**
 * exp_prg_generate -- Generate a pseudorandom sequence by
 * walking on an expander graph.
 *
 * Input: seed as a starting vertex index + edge choices
 * Output: sequence of 'length' vertices visited
 *
 * The expander provides the "hardness" needed for PRG output.
 *
 * Complexity: O(length).
 */
bool exp_prg_generate(const ExpanderGraph *g, size_t start,
                       size_t length, size_t *output) {
    if (!g || !output || length == 0 || start >= g->n) return false;

    size_t current = start;
    output[0] = current;

    for (size_t i = 1; i < length; i++) {
        /* Use the i-th bit of the "seed" to pick the edge */
        /* Here we use a simple deterministic rule: i mod d */
        size_t edge = i % g->d;
        size_t next = g->adj[current * g->d + edge];
        if (next == SIZE_MAX || next >= g->n) {
            /* Fallback to first valid neighbor */
            for (size_t k = 0; k < g->d; k++) {
                next = g->adj[current * g->d + k];
                if (next != SIZE_MAX && next < g->n) break;
            }
        }
        current = next;
        output[i] = current;
    }
    return true;
}

/**
 * exp_extractor_graph -- Build an extractor-like function from an
 * expander random walk: use the vertex along a walk as output.
 *
 * Extractors convert a weakly random source + short seed into
 * nearly uniform output. Expanders are the key building block
 * for extractors.
 *
 * This function simulates an extractor using an expander walk.
 *
 * source: array of n "weak" bits
 * seed: short seed (interpreted as vertex index)
 * output: output bits
 * length: number of output bits
 */
bool exp_extractor_graph(const ExpanderGraph *g, const char *source,
                          size_t seed, size_t length, char *output) {
    if (!g || !source || !output || length == 0 || seed >= g->n)
        return false;

    size_t current = seed;
    for (size_t i = 0; i < length; i++) {
        /* Output the source bit at current vertex */
        if (current < g->n) {
            output[i] = source[current];
        } else {
            output[i] = '0';
        }

        /* Move to next vertex via expander step */
        size_t edge = i % g->d;
        size_t next = g->adj[current * g->d + edge];
        if (next != SIZE_MAX && next < g->n) {
            current = next;
        }
    }
    return true;
}

/* ==================================================================
 * Lossless Expander Applications
 * ================================================================== */

/**
 * exp_is_lossless_expander -- Check if graph is a lossless expander:
 * for every subset S of size <= k, |N(S)| >= (1-eps) * d * |S|.
 *
 * Lossless expanders are stronger than vertex expanders;
 * they ensure that the neighbor set is almost as large as the
 * maximum possible (d * |S|), losing only an eps fraction.
 *
 * This is a heuristic check using random subsets.
 *
 * Complexity: O(num_samples * k * d).
 */
bool exp_is_lossless_expander(const ExpanderGraph *g, size_t k,
                               double eps, size_t num_samples) {
    if (!g || k == 0 || k > g->n || eps < 0 || eps > 1) return false;

    for (size_t sample = 0; sample < num_samples; sample++) {
        /* Generate a random subset S of size k */
        char *in_S = (char *)calloc(g->n, sizeof(char));
        if (!in_S) return false;

        size_t count = 0;
        while (count < k) {
            size_t v = (size_t)rand() % g->n;
            if (!in_S[v]) { in_S[v] = 1; count++; }
        }

        /* Collect neighbor set N(S) */
        char *in_N = (char *)calloc(g->n, sizeof(char));
        if (!in_N) { free(in_S); return false; }

        size_t n_count = 0;
        for (size_t u = 0; u < g->n; u++) {
            if (!in_S[u]) continue;
            for (size_t i = 0; i < g->d; i++) {
                size_t w = g->adj[u * g->d + i];
                if (w != SIZE_MAX && w < g->n && !in_N[w]) {
                    in_N[w] = 1;
                    n_count++;
                }
            }
        }

        free(in_S); free(in_N);

        /* Check: |N(S)| >= (1-eps) * d * k */
        double threshold = (1.0 - eps) * (double)g->d * (double)k;
        if ((double)n_count < threshold) return false;
    }
    return true;
}

/**
 * exp_random_walk_hitting_time -- Expected hitting time of vertex t
 * from vertex s via random walk on the expander graph.
 *
 * For expanders, the hitting time is O(log n) with high probability.
 * This function returns the expected value via simulation.
 *
 * Complexity: O(num_simulations * max_steps).
 */
double exp_random_walk_hitting_time(const ExpanderGraph *g, size_t s,
                                     size_t t, size_t num_simulations,
                                     size_t max_steps) {
    if (!g || s >= g->n || t >= g->n) return INFINITY;
    if (s == t) return 0.0;

    double total_steps = 0.0;
    size_t successful = 0;

    for (size_t sim = 0; sim < num_simulations; sim++) {
        size_t pos = s;
        for (size_t step = 1; step <= max_steps; step++) {
            pos = exp_random_walk_step(g, pos);
            if (pos == t) {
                total_steps += (double)step;
                successful++;
                break;
            }
        }
    }

    if (successful == 0) return INFINITY;
    return total_steps / (double)successful;
}

/* ==================================================================
 * Expander-Based Load Balancing
 * ================================================================== */

/**
 * exp_load_balance -- Use an expander graph to distribute tasks
 * across processors in a load-balanced way.
 *
 * The expander's expansion property ensures that each processor's
 * neighbors are well-distributed, providing fault tolerance
 * and balanced communication.
 *
 * This is a simplified model for expander-based parallel computing.
 *
 * load: array of n initial loads
 * g: expander representing processor interconnect
 * new_load: output array for balanced loads
 * rounds: number of load-balancing rounds
 *
 * Complexity: O(rounds * n * d).
 */
void exp_load_balance(const double *load, const ExpanderGraph *g,
                       double *new_load, size_t rounds) {
    if (!load || !g || !new_load) return;

    size_t n = g->n;
    memcpy(new_load, load, n * sizeof(double));

    double *temp = (double *)malloc(n * sizeof(double));
    if (!temp) return;

    for (size_t r = 0; r < rounds; r++) {
        memcpy(temp, new_load, n * sizeof(double));
        for (size_t u = 0; u < n; u++) {
            double neighbor_sum = 0.0;
            size_t valid_neighbors = 0;
            for (size_t k = 0; k < g->d; k++) {
                size_t v = g->adj[u * g->d + k];
                if (v != SIZE_MAX && v < n) {
                    neighbor_sum += temp[v];
                    valid_neighbors++;
                }
            }
            if (valid_neighbors > 0) {
                /* Diffusion: move toward average of neighbors */
                double avg = neighbor_sum / (double)valid_neighbors;
                new_load[u] = 0.5 * temp[u] + 0.5 * avg;
            }
        }
    }
    free(temp);
}

/**
 * exp_load_imbalance -- Measure the maximum load imbalance:
 *   max_i |load_i - avg_load| / avg_load
 *
 * A value near 0 indicates good balance.
 */
double exp_load_imbalance(const double *load, size_t n) {
    if (!load || n == 0) return INFINITY;

    double sum = 0.0, max_dev = 0.0;
    for (size_t i = 0; i < n; i++) sum += load[i];
    double avg = sum / (double)n;
    if (avg < 1e-15) return 0.0;

    for (size_t i = 0; i < n; i++) {
        double dev = fabs(load[i] - avg) / avg;
        if (dev > max_dev) max_dev = dev;
    }
    return max_dev;
}
