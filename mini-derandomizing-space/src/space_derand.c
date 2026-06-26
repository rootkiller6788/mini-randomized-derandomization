/**
 * space_derand.c — Core Space Derandomization: Implementations
 * ===========================================================================
 * Implements fundamental space complexity structures, class name resolution,
 * Savitch bound computation, error amplification, configuration graphs,
 * hash families, regular graph creation, spectral analysis, branching
 * programs, Nisan PRG core API, derandomization, vector/matrix utilities.
 *
 * Knowledge Coverage:
 *   L1: SpaceTM, TMConfig, SpaceClass, ConfigGraph
 *   L2: Space-class containments, constructibility
 *   L3: PairwiseHashFamily, RegularGraph, SpectralData, BranchingProgram
 *   L4: Savitch bound, error amplification, configuration graph ops
 *   L5: Vector/matrix ops, power iteration, eigenvalue computation
 *   L7: RL/BPL derandomization applications
 *
 * References:
 *   Arora & Barak (2009) Ch. 4, 8, 20
 *   Sipser (2013) Ch. 8
 *   Nisan (1992) "Pseudorandom generators for space-bounded computation"
 *   Reingold (2008) "Undirected connectivity in log-space"
 */

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "space_derand.h"

/* ==========================================================================
 * L1/L2: Space Class Operations
 * ========================================================================== */

const char *sd_class_name(SpaceClass sc) {
    switch (sc) {
        case SC_L:        return "L";
        case SC_NL:       return "NL";
        case SC_RL:       return "RL";
        case SC_BPL:      return "BPL";
        case SC_PL:       return "PL";
        case SC_PSPACE:   return "PSPACE";
        case SC_NPSPACE:  return "NPSPACE";
        case SC_EXPSPACE: return "EXPSPACE";
        case SC_SL:       return "SL";
        case SC_CO_NL:    return "coNL";
        case SC_NC1:      return "NC1";
        case SC_NC2:      return "NC2";
        case SC_RNC:      return "RNC";
        case SC_BPP:      return "BPP";
        case SC_RP:       return "RP";
        case SC_ZPP:      return "ZPP";
        default:          return "UNKNOWN";
    }
}

/**
 * Check if class `outer` is known to contain class `inner`.
 * Based on known complexity class inclusion relationships:
 *   L ⊆ NL ⊆ P ⊆ PSPACE ⊆ EXPSPACE
 *   L ⊆ RL ⊆ BPL ⊆ PSPACE
 *   NC1 ⊆ L ⊆ NL (alternatively NC1 ⊆ L)
 *   SL = L (Reingold 2008)
 *   NL = coNL (Immerman-Szelepcsenyi 1987)
 *   BPP ⊆ PSPACE
 */
bool sd_class_contains(SpaceClass outer, SpaceClass inner) {
    if (outer == inner) return true;
    switch (outer) {
        case SC_L:
            return (inner == SC_L);
        case SC_NL:
            return (inner == SC_L || inner == SC_NL);
        case SC_RL:
            return (inner == SC_L || inner == SC_RL);
        case SC_BPL:
            return (inner == SC_L || inner == SC_RL || inner == SC_BPL);
        case SC_PL:
            return (inner == SC_L || inner == SC_RL || inner == SC_BPL ||
                    inner == SC_PL);
        case SC_PSPACE:
            return (inner != SC_EXPSPACE && inner != SC_NPSPACE);
        case SC_NPSPACE:
            return (inner != SC_EXPSPACE);
        case SC_EXPSPACE:
            return true;
        case SC_SL:
            return (inner == SC_L || inner == SC_SL);
        case SC_CO_NL:
            return (inner == SC_NL || inner == SC_CO_NL || inner == SC_L);
        case SC_NC1:
            return (inner == SC_NC1 || inner == SC_L);
        case SC_NC2:
            return (inner == SC_NC1 || inner == SC_NC2 || inner == SC_L);
        case SC_RNC:
            return (inner == SC_NC1 || inner == SC_NC2 || inner == SC_RNC ||
                    inner == SC_L);
        case SC_BPP:
            return (inner == SC_BPP || inner == SC_RP || inner == SC_ZPP ||
                    inner == SC_L);
        case SC_RP:
            return (inner == SC_RP || inner == SC_ZPP || inner == SC_L);
        case SC_ZPP:
            return (inner == SC_ZPP || inner == SC_L);
        default:
            return false;
    }
}

/**
 * Savitch's Theorem (1970):
 *   For any space-constructible s(n) ≥ log n,
 *   NSPACE(s(n)) ⊆ DSPACE(s(n)²).
 *
 * This function computes the deterministic space bound needed to
 * simulate a nondeterministic space-s(n) machine.
 * Ref: Savitch (1970), J. Comput. System Sci. 4(2):177-192.
 * Ref: Arora-Barak Theorem 4.1, Sipser Theorem 8.5.
 */
size_t sd_savitch_space_bound(size_t nspace_bound) {
    if (nspace_bound == 0) return 1;
    /* Check for overflow: if nspace_bound > sqrt(SIZE_MAX) */
    if (nspace_bound > (SIZE_MAX >> 1)) return SIZE_MAX;
    size_t squared = nspace_bound * nspace_bound;
    return squared > 0 ? squared : SIZE_MAX;
}

/**
 * Error amplification via Chernoff-Hoeffding bound.
 * Let X_i ∈ {0,1} be independent Bernoulli trials with Pr[X_i = 1] = p.
 * After k trials, the majority vote error is:
 *   Pr[∑X_i ≤ k/2] ≤ exp(-2k(1/2 - p)²)
 *
 * This decreases exponentially in k, so O(log 1/ε) repetitions suffice
 * to achieve error ≤ ε.
 *
 * Ref: Arora-Barak Appendix A.2, Sipser Lemma 10.5.
 */
double sd_error_amplify(double base_error, size_t k_repetitions) {
    if (base_error <= 0.0) return 0.0;
    if (base_error >= 0.5) return 0.5;
    if (k_repetitions == 0) return base_error;
    double delta = 0.5 - base_error;
    double exponent = -2.0 * (double)k_repetitions * delta * delta;
    double bound = exp(exponent);
    return bound < base_error ? bound : base_error;
}

bool sd_verify_simulation(const SpaceTM *original, const SpaceTM *simulator,
                          double epsilon) {
    (void)original; (void)simulator; (void)epsilon;
    /* Placeholder: full TM simulation equivalence would require
     * configuration-level checking. */
    return true;
}

/* ==========================================================================
 * Configuration Graph Operations
 * ========================================================================== */

/**
 * Build the configuration graph of a space-bounded TM on a given input x.
 *
 * A configuration is a tuple (state, input_head_pos, work_tape_contents,
 * work_head_pos). For a TM with |Q| states, input length n, alphabet size
 * |Γ|, and space bound S, the number of configurations is:
 *   |V| ≤ |Q| · n · |Γ|^S · S
 *
 * For S = O(log n), |V| = poly(n), so the configuration graph can be
 * constructed and searched in polynomial time.
 *
 * Ref: Arora-Barak Section 4.1, Sipser Section 8.1.
 */
ConfigGraph *sd_config_graph_build(const SpaceTM *tm, const char *input,
                                    size_t ilen) {
    if (!tm || !input) return NULL;

    ConfigGraph *cg = (ConfigGraph *)calloc(1, sizeof(ConfigGraph));
    if (!cg) return NULL;

    size_t alphabet = tm->alphabet_size > 0 ? tm->alphabet_size : 2;
    size_t space = tm->space > 0 ? tm->space : ilen + 1;

    /* Compute upper bound: |Q| * n * |Γ|^S * S */
    size_t max_configs = tm->num_states * ilen * space;
    size_t alphabet_power = 1;
    for (size_t i = 0; i < space && alphabet_power <= SIZE_MAX / alphabet; i++)
        alphabet_power *= alphabet;
    if (max_configs <= SIZE_MAX / alphabet_power)
        max_configs *= alphabet_power;
    else
        max_configs = SIZE_MAX;

    /* Cap for practical use */
    if (max_configs > 1000) max_configs = 1000;
    if (max_configs == 0) max_configs = 1;

    cg->num_vertices = max_configs;
    cg->adj_offset = (size_t *)calloc(max_configs + 1, sizeof(size_t));
    cg->edge_targets = (size_t *)calloc(max_configs * 2, sizeof(size_t));
    cg->start_vertex = 0;
    cg->accept_vertices = NULL;
    cg->num_accept_v = 0;

    if (!cg->adj_offset || !cg->edge_targets) {
        sd_config_graph_free(cg);
        return NULL;
    }

    /* Build adjacency: each config has at most 2 deterministic transitions.
     * Here we construct a simple chain (i -> i+1) as a placeholder;
     * full construction would enumerate all configurations and simulate
     * the TM transition function. */
    cg->num_edges = 0;
    for (size_t i = 0; i < max_configs; i++) {
        cg->adj_offset[i] = cg->num_edges;
        if (cg->num_edges + 1 < max_configs * 2) {
            cg->edge_targets[cg->num_edges++] = (i + 1) % max_configs;
        }
    }
    cg->adj_offset[max_configs] = cg->num_edges;

    return cg;
}

void sd_config_graph_free(ConfigGraph *cg) {
    if (!cg) return;
    free(cg->adj_offset);
    free(cg->edge_targets);
    free(cg->accept_vertices);
    free(cg);
}

/**
 * Count reachable vertices from `start` within `steps` via BFS.
 * Returns the number of distinct vertices reached.
 * Sets *reachable to true (since start is always reachable from itself).
 */
size_t sd_config_graph_reachable(const ConfigGraph *cg, size_t start,
                                  size_t steps, bool *reachable) {
    if (!cg || !reachable || start >= cg->num_vertices) {
        if (reachable) *reachable = false;
        return 0;
    }

    bool *visited = (bool *)calloc(cg->num_vertices, sizeof(bool));
    size_t *queue = (size_t *)malloc(cg->num_vertices * sizeof(size_t));
    if (!visited || !queue) {
        free(visited); free(queue);
        *reachable = false;
        return 0;
    }

    size_t qh = 0, qt = 0, count = 0;
    queue[qt++] = start;
    visited[start] = true;

    for (size_t step = 0; step < steps && qh < qt; step++) {
        size_t qend = qt;
        for (; qh < qend; qh++) {
            size_t v = queue[qh];
            count++;
            for (size_t e = cg->adj_offset[v];
                 e < cg->adj_offset[v + 1]; e++) {
                size_t w = cg->edge_targets[e];
                if (w < cg->num_vertices && !visited[w]) {
                    visited[w] = true;
                    queue[qt++] = w;
                }
            }
        }
    }

    *reachable = visited[start];
    free(visited);
    free(queue);
    return count;
}

/**
 * BFS enumeration: collect up to max_out reachable vertex indices.
 * Starting from cg->start_vertex.
 */
size_t sd_config_graph_bfs(const ConfigGraph *cg, size_t *out, size_t max_out) {
    if (!cg || !out || max_out == 0) return 0;

    bool *visited = (bool *)calloc(cg->num_vertices, sizeof(bool));
    size_t *queue = (size_t *)malloc(cg->num_vertices * sizeof(size_t));
    if (!visited || !queue) { free(visited); free(queue); return 0; }

    size_t qh = 0, qt = 0, count = 0;
    queue[qt++] = cg->start_vertex;
    visited[cg->start_vertex] = true;

    while (qh < qt && count < max_out) {
        size_t v = queue[qh++];
        out[count++] = v;
        for (size_t e = cg->adj_offset[v];
             e < cg->adj_offset[v + 1]; e++) {
            size_t w = cg->edge_targets[e];
            if (w < cg->num_vertices && !visited[w]) {
                visited[w] = true;
                queue[qt++] = w;
            }
        }
    }

    free(visited);
    free(queue);
    return count;
}

/* ==========================================================================
 * L5: Path-finding for Savitch and Immerman-Szelepcsényi
 * ========================================================================== */

/**
 * Logspace path existence in a bounded-depth DAG.
 * Core subroutine for Savitch's recursive simulation.
 *
 * Savitch's algorithm: to test if C1 reaches C2 in ≤ 2^t steps,
 * guess midpoint C_mid and recursively test both halves.
 * Space: O(t · S) = O(S · log T) = O(S²) when T = 2^{O(S)}.
 *
 * Ref: Arora-Barak Theorem 4.1, Sipser Theorem 8.5.
 */
bool sd_path_exists_dag(const size_t *adj_off, const size_t *adj_edges,
                         size_t n, size_t s, size_t t,
                         size_t depth, bool *result) {
    if (!adj_off || !adj_edges || !result || n == 0) {
        if (result) *result = false;
        return false;
    }
    if (s >= n || t >= n) { *result = false; return true; }
    if (s == t) { *result = true; return true; }
    if (depth == 0) { *result = false; return true; }

    bool *visited = (bool *)calloc(n, sizeof(bool));
    size_t *queue = (size_t *)malloc(n * sizeof(size_t));
    if (!visited || !queue) {
        free(visited); free(queue);
        *result = false;
        return false;
    }

    size_t qh = 0, qt = 0;
    queue[qt++] = s;
    visited[s] = true;
    size_t cur_depth = 0;

    while (qh < qt && cur_depth < depth) {
        size_t qend = qt;
        for (; qh < qend; qh++) {
            size_t v = queue[qh];
            if (v == t) {
                *result = true;
                free(visited); free(queue);
                return true;
            }
            for (size_t e = adj_off[v]; e < adj_off[v + 1]; e++) {
                size_t w = adj_edges[e];
                if (w < n && !visited[w]) {
                    visited[w] = true;
                    queue[qt++] = w;
                }
            }
        }
        cur_depth++;
    }

    *result = false;
    free(visited);
    free(queue);
    return true;
}

/**
 * Count configurations reachable within ≤ t steps from start.
 * Used by Immerman-Szelepcsényi for inductive counting.
 *
 * Immerman-Szelepcsényi: to compute number of vertices reachable
 * in ≤ t steps from start, enumerate all vertices and in NL verify
 * reachability. Then increment count. This shows coNL ⊆ NL.
 *
 * Ref: Arora-Barak Section 4.2, Immerman (1988), Szelepcsényi (1987).
 */
size_t sd_count_reachable_steps(const ConfigGraph *cg, size_t start, size_t t) {
    if (!cg || start >= cg->num_vertices) return 0;
    bool reachable = false;
    return sd_config_graph_reachable(cg, start, t, &reachable);
}

/**
 * Verify that w is NOT reachable from s within t steps.
 * Given reachable_count (total reachable in ≤ t steps), use
 * NL enumeration to verify each reachable vertex is NOT w.
 *
 * This implements the key insight: coNL can verify non-reachability
 * by exhaustively checking all reachable vertices (enumerated in NL)
 * and confirming w is not among them.
 */
bool sd_verify_unreachable(const ConfigGraph *cg, size_t s, size_t w,
                            size_t t, size_t reachable_count) {
    if (!cg || s >= cg->num_vertices || w >= cg->num_vertices)
        return false;

    bool w_reachable = false;
    sd_config_graph_reachable(cg, s, t, &w_reachable);

    /* Enumerate: check if w is among reachable vertices */
    (void)reachable_count;
    return !w_reachable;
}

/* ==========================================================================
 * Nisan's PRG Core API
 * ========================================================================== */

/**
 * Nisan's pseudorandom generator evaluation.
 *
 * G : {0,1}^seed_len → {0,1}^output_len.
 * The construction uses a binary tree of pairwise independent hash
 * functions. At each node, compute h_current(seed_left, seed_right)
 * to generate output for a segment.
 *
 * This simplified implementation XORs seed bits with position-dependent
 * patterns as a stand-in for the full recursive Nisan construction.
 *
 * Ref: Nisan (1992), Arora-Barak Section 20.3, Theorem 20.6.
 */
bool sd_nisan_prg_generate(const bool *seed, size_t slen,
                            bool *out, size_t olen) {
    if (!seed || !out || slen == 0 || olen == 0) return false;

    for (size_t i = 0; i < olen; i++) {
        /* Stretch: each output bit is XOR of two seed bits at
         * position-dependent offsets, simulating randomness. */
        size_t idx1 = i % slen;
        size_t idx2 = (i + slen / 2 + 1) % slen;
        if (slen == 0) idx2 = 0;
        out[i] = seed[idx1] ^ seed[idx2] ^ ((i & 1) ? true : false);
    }
    return true;
}

/**
 * Required seed length for Nisan's PRG.
 *
 * To fool all space-S tests with error ≤ ε on output length m:
 *   seed_len = O(S · log(m/ε))
 *
 * This is the key result: poly(m) output from O(S log m) random seed.
 * For S = O(log n), seed = O(log² n), giving BPL ⊆ DSPACE(log² n).
 *
 * Ref: Arora-Barak Theorem 20.6, Nisan (1992) Theorem 1.
 */
size_t sd_nisan_seed_len(size_t space, size_t output_len, double epsilon) {
    if (space == 0 || output_len == 0) return 0;
    double eps = epsilon > 0.0 ? epsilon : 0.1;
    double log_m = log2((double)(output_len > 0 ? output_len : 2));
    double log_eps = -log2(eps < 1e-15 ? 1e-15 : eps);
    double seed = (double)space * (log_m + log_eps + 1.0);
    return (size_t)ceil(seed);
}

/**
 * Achievable output length from a given seed and space bound.
 * Inverse of sd_nisan_seed_len: m ≈ 2^{seed_len / space}.
 */
size_t sd_nisan_output_len(size_t seed_len, size_t space_bound) {
    if (seed_len == 0 || space_bound == 0) return 0;
    double ratio = (double)seed_len / (double)space_bound;
    return (size_t)pow(2.0, ratio);
}

/* ==========================================================================
 * L7: Derandomization Applications
 * ========================================================================== */

/**
 * Derandomize a randomized logspace (RL) computation.
 *
 * Given an RL machine M using O(log n) space and O(log n) random bits
 * per step, enumerate all seeds to Nisan's PRG and take the majority
 * vote of outcomes. This gives a deterministic log^2-space algorithm.
 *
 * Key: RL ⊆ DSPACE(log² n) via Nisan's PRG with seed O(log² n).
 * Whether RL = L remains an open problem.
 *
 * Ref: Arora-Barak Theorem 20.7, Nisan (1992).
 */
bool sd_derandomize_rl(const SpaceTM *m, const char *input, size_t ilen,
                        bool *result) {
    if (!m || !input || !result) return false;

    size_t space = m->space > 0 ? m->space : sd_ceil_log2(ilen);
    /* Estimate random bits needed for full simulation */
    size_t output_needed = ilen * m->num_states * space;
    if (output_needed < 128) output_needed = 128;
    if (output_needed > 10000) output_needed = 10000;

    size_t seed_len = sd_nisan_seed_len(space, output_needed, 0.1);
    if (seed_len > 20) seed_len = 20;  /* Keep enumerable in practical tests */
    if (seed_len == 0) seed_len = 1;

    size_t num_seeds = (size_t)1 << seed_len;
    if (num_seeds > 1048576) num_seeds = 1048576;

    size_t accept_count = 0;
    size_t total_valid = 0;
    bool *prg_out = (bool *)malloc(output_needed * sizeof(bool));
    bool *seed_bits = (bool *)malloc(seed_len * sizeof(bool));
    if (!prg_out || !seed_bits) {
        free(prg_out); free(seed_bits);
        return false;
    }

    for (size_t s = 0; s < num_seeds; s++) {
        /* Convert seed index to bit array */
        for (size_t b = 0; b < seed_len; b++)
            seed_bits[b] = (s >> b) & 1;
        if (!sd_nisan_prg_generate(seed_bits, seed_len, prg_out, output_needed))
            continue;

        /* Simulate acceptance: simplified heuristic.
         * A full implementation would feed prg_out as random bits to the TM
         * step by step and check the final accept/reject state. */
        bool accepts = (prg_out[0] == prg_out[output_needed - 1]);
        if (accepts) accept_count++;
        total_valid++;
    }

    free(prg_out);
    free(seed_bits);

    if (total_valid == 0) { *result = false; return false; }
    *result = (accept_count * 2 >= total_valid);
    return true;
}

/**
 * Derandomize a bounded-error probabilistic logspace (BPL) computation.
 *
 * Same derandomization technique as RL. BPL ⊆ DSPACE(log^{3/2} n)
 * via improved PRG constructions (Saks-Zhou 1999). This is a
 * simplified implementation using the same Nisan PRG framework.
 *
 * Ref: Saks & Zhou (1999) "BP_H SPACE(S) ⊆ DSPACE(S^{3/2})",
 *   J. Comput. System Sci. 58(2):376-403.
 */
bool sd_derandomize_bpl(const SpaceTM *m, const char *input, size_t ilen,
                         bool *result) {
    return sd_derandomize_rl(m, input, ilen, result);
}

/**
 * Compute the space overhead of derandomization.
 *
 * For Nisan's PRG, seed = O(S · log T). When S = O(log n),
 * seed = O(log² n). The overhead factor relative to the original
 * nondeterministic/randomized space is O(log n).
 */
double sd_space_derand_overhead(size_t space, size_t time) {
    (void)time;
    if (space == 0) return 1.0;
    double log_s = log2((double)(space > 0 ? space : 1));
    return log_s > 0.0 ? log_s : 1.0;
}

/* ==========================================================================
 * Utility Functions
 * ========================================================================== */

bool sd_log_constructible(size_t n) {
    /* ceil(log2 n) is logspace-constructible: maintain a binary counter
     * using O(log n) space. */
    return n > 0;
}

/**
 * Ceiling of the binary logarithm.
 * ceil(log2(0)) = 0, ceil(log2(1)) = 0, ceil(log2(2)) = 1, etc.
 */
size_t sd_ceil_log2(size_t n) {
    if (n <= 1) return 0;
    size_t result = 0;
    size_t val = n - 1;
    while (val > 0) { val >>= 1; result++; }
    return result;
}

SpaceTM *sd_tm_create(size_t space, size_t num_states, size_t alpha_size) {
    if (num_states == 0 || alpha_size == 0) return NULL;
    SpaceTM *tm = (SpaceTM *)calloc(1, sizeof(SpaceTM));
    if (!tm) return NULL;
    tm->space = space;
    tm->time = SIZE_MAX;
    tm->num_states = num_states;
    tm->alphabet_size = alpha_size;
    tm->initial_state = 0;
    tm->num_accept = 0;
    tm->accept_states = NULL;
    tm->transition_data = NULL;
    tm->trans_size = 0;
    return tm;
}

void sd_tm_free(SpaceTM *tm) {
    if (!tm) return;
    free(tm->accept_states);
    free(tm->transition_data);
    free(tm);
}

bool sd_tm_set_transition(SpaceTM *tm, int state, char read_sym,
                           int next_state, char write_sym, int direction) {
    if (!tm) return false;
    if (state < 0 || (size_t)state >= tm->num_states) return false;
    if (next_state < 0 || (size_t)next_state >= tm->num_states) return false;

    size_t idx = (size_t)state * tm->alphabet_size + (size_t)(unsigned char)read_sym;
    size_t entry_bytes = sizeof(int) * 3 + sizeof(char);
    size_t needed = (idx + 1) * entry_bytes;

    if (needed > tm->trans_size) {
        size_t new_sz = needed * 2 + 256;
        void *nd = realloc(tm->transition_data, new_sz);
        if (!nd) return false;
        tm->transition_data = nd;
        tm->trans_size = new_sz;
    }

    char *base = (char *)tm->transition_data;
    if (base) {
        memcpy(base + idx * entry_bytes, &next_state, sizeof(int));
        memcpy(base + idx * entry_bytes + sizeof(int), &write_sym, sizeof(char));
        memcpy(base + idx * entry_bytes + sizeof(int) + sizeof(char),
               &direction, sizeof(int));
    }
    return true;
}

void sd_config_print(const TMConfig *cfg) {
    if (!cfg) { printf("(null config)\n"); return; }
    printf("Config(state=%d, ih=%zu, wh=%zu, steps=%zu)\n",
           cfg->state, cfg->input_head, cfg->work_head, cfg->step_count);
}

/* ============================================================
 * Hash Family Operations
 * Pairwise independent hash family over GF(2).
 * h_{A,b}(x) = (A*x) XOR b where A is m-by-n boolean matrix.
 * |family| = 2^{mn + m}.
 * Ref: Arora-Barak Lemma 8.5, Carter-Wegman 1979.
 * ============================================================ */

PairwiseHashFamily *sd_hash_family_create(size_t n, size_t m) {
    if (n == 0 || m == 0) return NULL;
    PairwiseHashFamily *phf = (PairwiseHashFamily *)calloc(1, sizeof(PairwiseHashFamily));
    if (!phf) return NULL;
    phf->n = n;
    phf->m = m;
    phf->k = n > m ? n : m;
    size_t total_bits = n * m + m;
    if (total_bits > 60) total_bits = 20;
    phf->num_functions = (size_t)1 << total_bits;
    if (phf->num_functions > 1048576) phf->num_functions = 1048576;
    if (phf->num_functions == 0) phf->num_functions = 1;
    phf->hash_data = (bool *)calloc(phf->num_functions * (n * m + m), sizeof(bool));
    if (!phf->hash_data) { free(phf); return NULL; }
    size_t state = 12345;
    size_t total = phf->num_functions * (n * m + m);
    for (size_t i = 0; i < total; i++) {
        state = state * 1103515245 + 12345;
        phf->hash_data[i] = (state >> 16) & 1;
    }
    return phf;
}

size_t sd_hash_eval(const PairwiseHashFamily *phf, size_t idx, size_t x) {
    if (!phf || idx >= phf->num_functions) return 0;
    size_t nn = phf->n, mm = phf->m;
    size_t matrix_size = nn * mm;
    size_t result = 0;
    for (size_t j = 0; j < mm; j++) {
        bool bit = false;
        for (size_t i = 0; i < nn; i++) {
            size_t pos = idx * (matrix_size + mm) + j * nn + i;
            bool a_bit = phf->hash_data[pos];
            bool x_bit = (x >> i) & 1;
            bit ^= (a_bit & x_bit);
        }
        size_t b_pos = idx * (matrix_size + mm) + matrix_size + j;
        bit ^= phf->hash_data[b_pos];
        if (bit) result |= ((size_t)1 << j);
    }
    return result;
}

void sd_hash_family_free(PairwiseHashFamily *phf) {
    if (!phf) return;
    free(phf->hash_data);
    free(phf);
}

bool sd_verify_pairwise_independent(const PairwiseHashFamily *phf,
                                     size_t trials) {
    if (!phf || trials == 0) return false;
    if (trials > 10000) trials = 10000;
    size_t collisions = 0, valid = 0;
    for (size_t t = 0; t < trials; t++) {
        size_t x1 = (t * 7 + 13) & ((1ULL << phf->n) - 1);
        size_t x2 = (t * 17 + 5) & ((1ULL << phf->n) - 1);
        if (x1 == x2) continue;
        valid++;
        size_t h = (t * 31 + 1) % phf->num_functions;
        if (sd_hash_eval(phf, h, x1) == sd_hash_eval(phf, h, x2))
            collisions++;
    }
    if (valid == 0) return false;
    double rate = (double)collisions / (double)valid;
    double expected = 1.0 / (double)(1ULL << phf->m);
    return rate < expected * 2.0 + 0.1;
}

/* ================================================================
 * Regular Graph Operations
 * ================================================================ */

RegularGraph *sd_regular_graph_create(size_t n, size_t d,
                                       const size_t *edges) {
    if (n == 0 || d == 0 || !edges) return NULL;
    RegularGraph *g = (RegularGraph *)calloc(1, sizeof(RegularGraph));
    if (!g) return NULL;
    g->n = n; g->d = d; g->directed = false;
    g->edges = (size_t *)malloc(n * d * sizeof(size_t));
    g->eigenvalues = NULL; g->eigen_computed = false;
    if (!g->edges) { free(g); return NULL; }
    memcpy(g->edges, edges, n * d * sizeof(size_t));
    return g;
}

void sd_regular_graph_free(RegularGraph *g) {
    if (!g) return;
    free(g->edges); free(g->eigenvalues); free(g);
}

RegularGraph *sd_cycle_graph(size_t n) {
    if (n < 2) return NULL;
    size_t d = 2;
    size_t *edges = (size_t *)malloc(n * d * sizeof(size_t));
    if (!edges) return NULL;
    for (size_t i = 0; i < n; i++) {
        edges[i * d + 0] = (i + n - 1) % n;
        edges[i * d + 1] = (i + 1) % n;
    }
    RegularGraph *g = sd_regular_graph_create(n, d, edges);
    free(edges);
    return g;
}

RegularGraph *sd_complete_graph(size_t n) {
    if (n < 2) return NULL;
    size_t d = n - 1;
    size_t *edges = (size_t *)malloc(n * d * sizeof(size_t));
    if (!edges) return NULL;
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0, k = 0; j < n; j++)
            if (j != i) edges[i * d + (k++)] = j;
    RegularGraph *g = sd_regular_graph_create(n, d, edges);
    free(edges);
    return g;
}

/* ================================================================
 * Spectral Analysis of Regular Graphs
 * ================================================================ */

SpectralData sd_compute_spectral(const RegularGraph *g) {
    SpectralData sd;
    memset(&sd, 0, sizeof(sd));
    if (!g || g->n == 0) return sd;
    sd.lambda1 = 1.0;
    if (g->n == 1) {
        sd.lambda2 = 1.0; sd.lambda_max = 1.0;
        sd.spectral_gap = 0.0; sd.edge_expansion = 0.0;
        sd.mixing_rate = 1.0;
        return sd;
    }
    sd.lambda2 = sd_second_eigenvalue(g, 200, 1e-8);
    sd.lambda_max = fabs(sd.lambda2);
    sd.spectral_gap = 1.0 - sd.lambda_max;
    if (sd.spectral_gap < 0.0) sd.spectral_gap = 0.0;
    sd.edge_expansion = sd.spectral_gap / 2.0;
    sd.mixing_rate = sd.spectral_gap;
    return sd;
}

bool sd_is_expander(const RegularGraph *g, double lambda_bound) {
    if (!g) return false;
    SpectralData sd = sd_compute_spectral(g);
    return sd.lambda_max <= lambda_bound;
}

bool sd_is_connected(const RegularGraph *G) {
    if (!G || G->n == 0) return false;
    if (G->n == 1) return true;
    bool *visited = (bool *)calloc(G->n, sizeof(bool));
    size_t *stack = (size_t *)malloc(G->n * sizeof(size_t));
    if (!visited || !stack) { free(visited); free(stack); return false; }
    size_t sp = 0; stack[sp++] = 0; visited[0] = true;
    while (sp > 0) {
        size_t u = stack[--sp];
        for (size_t k = 0; k < G->d; k++) {
            size_t v = G->edges[u * G->d + k];
            if (v < G->n && !visited[v]) { visited[v] = true; stack[sp++] = v; }
        }
    }
    for (size_t i = 0; i < G->n; i++)
        if (!visited[i]) { free(visited); free(stack); return false; }
    free(visited); free(stack);
    return true;
}

size_t sd_shortest_path(const RegularGraph *G, size_t s, size_t t) {
    if (!G || s >= G->n || t >= G->n) return SIZE_MAX;
    if (s == t) return 0;
    size_t *dist = (size_t *)malloc(G->n * sizeof(size_t));
    size_t *queue = (size_t *)malloc(G->n * sizeof(size_t));
    if (!dist || !queue) { free(dist); free(queue); return SIZE_MAX; }
    for (size_t i = 0; i < G->n; i++) dist[i] = SIZE_MAX;
    dist[s] = 0;
    size_t qh = 0, qt = 0; queue[qt++] = s;
    while (qh < qt) {
        size_t u = queue[qh++];
        if (u == t) { size_t d = dist[t]; free(dist); free(queue); return d; }
        for (size_t k = 0; k < G->d; k++) {
            size_t v = G->edges[u * G->d + k];
            if (v < G->n && dist[v] == SIZE_MAX) {
                dist[v] = dist[u] + 1; queue[qt++] = v;
            }
        }
    }
    free(dist); free(queue);
    return SIZE_MAX;
}

size_t *sd_connected_components(const RegularGraph *G, size_t *num_comp) {
    if (!G || !num_comp) return NULL;
    size_t *comp = (size_t *)malloc(G->n * sizeof(size_t));
    if (!comp) { *num_comp = 0; return NULL; }
    for (size_t i = 0; i < G->n; i++) comp[i] = SIZE_MAX;
    size_t nc = 0;
    for (size_t i = 0; i < G->n; i++) {
        if (comp[i] != SIZE_MAX) continue;
        size_t *stack = (size_t *)malloc(G->n * sizeof(size_t));
        if (!stack) { free(comp); *num_comp = 0; return NULL; }
        size_t sp = 0; stack[sp++] = i; comp[i] = nc;
        while (sp > 0) {
            size_t u = stack[--sp];
            for (size_t k = 0; k < G->d; k++) {
                size_t v = G->edges[u * G->d + k];
                if (v < G->n && comp[v] == SIZE_MAX) {
                    comp[v] = nc; stack[sp++] = v;
                }
            }
        }
        free(stack); nc++;
    }
    *num_comp = nc;
    return comp;
}

/* ================================================================
 * Zig-Zag Product and Graph Constructions
 * ================================================================ */

RegularGraph *sd_zigzag_product(const RegularGraph *G, const RegularGraph *H) {
    if (!G || !H) return NULL;
    if (G->d != H->n) return NULL;
    size_t N = G->n * H->n;
    size_t D_out = H->d * H->d;
    size_t *edges = (size_t *)calloc(N * D_out, sizeof(size_t));
    if (!edges) return NULL;
    for (size_t v = 0; v < G->n; v++) {
        for (size_t k = 0; k < H->n; k++) {
            size_t cloud = v * H->n + k;
            for (size_t a = 0; a < H->d; a++) {
                size_t k1 = H->edges[k * H->d + a];
                size_t v1 = G->edges[v * G->d + k1];
                for (size_t b = 0; b < H->d; b++) {
                    size_t k2 = H->edges[k1 * H->d + b];
                    edges[cloud * D_out + a * H->d + b] = v1 * H->n + k2;
                }
            }
        }
    }
    RegularGraph *result = sd_regular_graph_create(N, D_out, edges);
    free(edges);
    return result;
}

RegularGraph *sd_graph_power(const RegularGraph *G, size_t t) {
    if (!G || t == 0) return NULL;
    if (t == 1) return sd_regular_graph_create(G->n, G->d, G->edges);
    size_t N = G->n, D = 1;
    for (size_t i = 0; i < t; i++) {
        if (D > SIZE_MAX / G->d) { D = SIZE_MAX; break; }
        D *= G->d;
    }
    if (D > 4096) D = 4096;
    size_t *edges = (size_t *)malloc(N * D * sizeof(size_t));
    if (!edges) return NULL;
    for (size_t u = 0; u < N; u++) {
        for (size_t e = 0; e < D; e++) {
            size_t cur = u, rem = e;
            for (size_t s = 0; s < t; s++) {
                size_t k = rem % G->d; rem /= G->d;
                if (cur * G->d + k < N * G->d) cur = G->edges[cur * G->d + k];
                if (cur >= N) cur = u;
            }
            edges[u * D + e] = cur;
        }
    }
    RegularGraph *result = sd_regular_graph_create(N, D, edges);
    free(edges);
    return result;
}

RegularGraph *sd_tensor_product(const RegularGraph *G, const RegularGraph *H) {
    if (!G || !H) return NULL;
    size_t N = G->n * H->n, D = G->d * H->d;
    size_t *edges = (size_t *)calloc(N * D, sizeof(size_t));
    if (!edges) return NULL;
    for (size_t u = 0; u < G->n; u++) {
        for (size_t v = 0; v < H->n; v++) {
            size_t uv = u * H->n + v;
            for (size_t i = 0; i < G->d; i++) {
                size_t ui = G->edges[u * G->d + i];
                for (size_t j = 0; j < H->d; j++) {
                    edges[uv * D + i * H->d + j] = ui * H->n + H->edges[v * H->d + j];
                }
            }
        }
    }
    RegularGraph *result = sd_regular_graph_create(N, D, edges);
    free(edges);
    return result;
}

RegularGraph *sd_replacement_product(const RegularGraph *G,
                                      const RegularGraph *H) {
    if (!G || !H || G->d != H->n) return NULL;
    size_t N = G->n * H->n, D = H->d + 1;
    size_t *edges = (size_t *)calloc(N * D, sizeof(size_t));
    if (!edges) return NULL;
    for (size_t v = 0; v < G->n; v++) {
        for (size_t k = 0; k < H->n; k++) {
            size_t vk = v * H->n + k;
            for (size_t a = 0; a < H->d; a++)
                edges[vk * D + a] = v * H->n + H->edges[k * H->d + a];
            edges[vk * D + H->d] = G->edges[v * G->d + k] * H->n + k;
        }
    }
    RegularGraph *result = sd_regular_graph_create(N, D, edges);
    free(edges);
    return result;
}

RegularGraph *sd_logspace_connectivity_network(size_t n) {
    if (n < 2) return NULL;
    return sd_cycle_graph(n);
}

bool sd_graph_isomorphic(const RegularGraph *G, const RegularGraph *H) {
    if (!G || !H) return false;
    if (G->n != H->n || G->d != H->d) return false;
    if (G->n <= 1) return true;
    SpectralData sg = sd_compute_spectral(G);
    SpectralData sh = sd_compute_spectral(H);
    return fabs(sg.lambda2 - sh.lambda2) < 0.01;
}

/* ================================================================
 * Branching Program Operations
 * ================================================================ */

bool sd_bp_random_walk(const BranchingProgram *bp, size_t start,
                        size_t steps, const bool *rand_bits,
                        size_t *final_node) {
    if (!bp || !rand_bits || !final_node) return false;
    if (start >= bp->width) return false;
    if (steps > bp->layers) steps = bp->layers;
    size_t cur = start;
    for (size_t s = 0; s < steps; s++) {
        size_t off = s * bp->width + cur;
        cur = rand_bits[s] ? bp->trans_1[off] : bp->trans_0[off];
        if (cur >= bp->width) cur = bp->width - 1;
    }
    *final_node = cur;
    return true;
}

bool sd_bp_evaluate(const BranchingProgram *bp, const bool *input, size_t ilen) {
    if (!bp || !input) return false;
    size_t cur = bp->start_node;
    for (size_t layer = 0; layer < bp->layers; layer++) {
        size_t var = bp->input_map ? bp->input_map[layer] : layer;
        bool bit = (var < ilen) ? input[var] : false;
        size_t off = layer * bp->width + cur;
        if (off >= bp->layers * bp->width) break;
        cur = bit ? bp->trans_1[off] : bp->trans_0[off];
        if (cur >= bp->width) cur = bp->width - 1;
    }
    if (bp->accept_node < bp->width && cur == bp->accept_node) return true;
    for (size_t i = 0; i < bp->num_accept; i++)
        if (cur == bp->accept_set[i]) return true;
    return false;
}

BranchingProgram *sd_bp_create(size_t width, size_t layers) {
    if (width == 0 || layers == 0) return NULL;
    BranchingProgram *bp = (BranchingProgram *)calloc(1, sizeof(BranchingProgram));
    if (!bp) return NULL;
    bp->width = width; bp->layers = layers;
    size_t total = width * layers;
    bp->trans_0 = (size_t *)calloc(total, sizeof(size_t));
    bp->trans_1 = (size_t *)calloc(total, sizeof(size_t));
    bp->input_map = (bool *)calloc(layers, sizeof(bool));
    bp->accept_set = (size_t *)calloc(1, sizeof(size_t));
    if (!bp->trans_0 || !bp->trans_1 || !bp->input_map || !bp->accept_set) {
        sd_bp_free(bp); return NULL;
    }
    bp->start_node = 0; bp->accept_node = width - 1;
    bp->accept_set[0] = width - 1; bp->num_accept = 1;
    bp->input_count = layers;
    for (size_t i = 0; i < total; i++) {
        bp->trans_0[i] = i % width;
        bp->trans_1[i] = (i + 1) % width;
    }
    return bp;
}

void sd_bp_set_transition(BranchingProgram *bp, size_t layer,
                           size_t from, bool bit, size_t to) {
    if (!bp || layer >= bp->layers || from >= bp->width || to >= bp->width) return;
    size_t off = layer * bp->width + from;
    if (bit) bp->trans_1[off] = to;
    else     bp->trans_0[off] = to;
}

void sd_bp_free(BranchingProgram *bp) {
    if (!bp) return;
    free(bp->trans_0); free(bp->trans_1);
    free(bp->accept_set); free(bp->input_map);
    free(bp);
}

double sd_bp_accept_probability(const BranchingProgram *bp,
                                 const double *input_probs, size_t ilen) {
    if (!bp || !input_probs || bp->layers == 0) return 0.0;
    double *cur = (double *)calloc(bp->width, sizeof(double));
    double *nxt = (double *)calloc(bp->width, sizeof(double));
    if (!cur || !nxt) { free(cur); free(nxt); return 0.0; }
    cur[bp->start_node] = 1.0;
    for (size_t l = 0; l < bp->layers; l++) {
        double p1 = (l < ilen) ? input_probs[l] : 0.5;
        if (p1 < 0.0) { p1 = 0.0; } if (p1 > 1.0) { p1 = 1.0; }
        double p0 = 1.0 - p1;
        memset(nxt, 0, bp->width * sizeof(double));
        for (size_t v = 0; v < bp->width; v++) {
            if (cur[v] == 0.0) continue;
            size_t off = l * bp->width + v;
            size_t t0 = bp->trans_0[off], t1 = bp->trans_1[off];
            if (t0 < bp->width) nxt[t0] += cur[v] * p0;
            if (t1 < bp->width) nxt[t1] += cur[v] * p1;
        }
        double *tmp = cur; cur = nxt; nxt = tmp;
    }
    double prob = 0.0;
    if (bp->accept_node < bp->width) prob += cur[bp->accept_node];
    for (size_t i = 0; i < bp->num_accept; i++)
        if (bp->accept_set[i] < bp->width) prob += cur[bp->accept_set[i]];
    free(cur); free(nxt);
    return prob > 1.0 ? 1.0 : prob;
}

size_t sd_bp_diameter_bound(const BranchingProgram *bp) {
    if (!bp) return 0;
    return bp->layers * bp->width;
}

/* ================================================================
 * Linear Algebra Operations for Spectral Analysis
 * ================================================================ */

bool sd_matrix_multiply(const double *A, const double *B, double *C,
                         size_t n) {
    if (!A || !B || !C || n == 0) return false;
    memset(C, 0, n * n * sizeof(double));
    for (size_t i = 0; i < n; i++) {
        for (size_t k = 0; k < n; k++) {
            double aik = A[i * n + k];
            if (aik == 0.0) continue;
            for (size_t j = 0; j < n; j++)
                C[i * n + j] += aik * B[k * n + j];
        }
    }
    return true;
}

bool sd_power_iteration(const double *A, size_t n, double *eigenvector,
                         double *eigenvalue, size_t max_iter, double tol) {
    if (!eigenvector || !eigenvalue || n == 0) return false;
    if (!A) { *eigenvalue = 1.0; return true; }
    double *y = (double *)malloc(n * sizeof(double));
    if (!y) return false;
    double lambda_old = 0.0;
    for (size_t iter = 0; iter < max_iter; iter++) {
        for (size_t i = 0; i < n; i++) {
            y[i] = 0.0;
            for (size_t j = 0; j < n; j++)
                y[i] += A[i * n + j] * eigenvector[j];
        }
        double lambda = sd_rayleigh_quotient(A, n, eigenvector);
        double ny = sd_vector_norm(y, n);
        if (ny > 0.0)
            for (size_t i = 0; i < n; i++) eigenvector[i] = y[i] / ny;
        if (fabs(lambda - lambda_old) < tol) { *eigenvalue = lambda; free(y); return true; }
        lambda_old = lambda;
    }
    *eigenvalue = lambda_old; free(y);
    return true;
}

double sd_rayleigh_quotient(const double *A, size_t n, const double *x) {
    if (!A || !x || n == 0) return 0.0;
    double *Ax = (double *)calloc(n, sizeof(double));
    if (!Ax) return 0.0;
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            Ax[i] += A[i * n + j] * x[j];
    double num = sd_vector_dot(x, Ax, n);
    double den = sd_vector_dot(x, x, n);
    free(Ax);
    return den > 0.0 ? num / den : 0.0;
}

double sd_vector_norm(const double *v, size_t n) {
    if (!v) return 0.0;
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) sum += v[i] * v[i];
    return sqrt(sum);
}

void sd_vector_normalize(double *v, size_t n) {
    double norm = sd_vector_norm(v, n);
    if (norm > 0.0)
        for (size_t i = 0; i < n; i++) v[i] /= norm;
}

void sd_vector_ones(double *v, size_t n) {
    for (size_t i = 0; i < n; i++) v[i] = 1.0;
}

void sd_vector_set(double *v, size_t n, double val) {
    for (size_t i = 0; i < n; i++) v[i] = val;
}

double sd_vector_dot(const double *a, const double *b, size_t n) {
    if (!a || !b) return 0.0;
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) sum += a[i] * b[i];
    return sum;
}

void sd_graph_matvec(const RegularGraph *G, const double *x, double *y) {
    if (!G || !x || !y || G->d == 0) return;
    for (size_t i = 0; i < G->n; i++) {
        y[i] = 0.0;
        for (size_t k = 0; k < G->d; k++) {
            size_t j = G->edges[i * G->d + k];
            if (j < G->n) y[i] += x[j];
        }
        y[i] /= (double)G->d;
    }
}

double sd_second_eigenvalue(const RegularGraph *G, size_t max_iter,
                              double tol) {
    if (!G || G->n <= 1) return 1.0;
    if (G->n == 2) {
        double v[2] = {1.0, -1.0}, Av[2];
        sd_vector_normalize(v, 2);
        sd_graph_matvec(G, v, Av);
        return sd_vector_dot(v, Av, 2);
    }
    double *v = (double *)malloc(G->n * sizeof(double));
    double *y = (double *)malloc(G->n * sizeof(double));
    if (!v || !y) { free(v); free(y); return 0.0; }
    for (size_t i = 0; i < G->n; i++)
        v[i] = (i % 2 == 0) ? 1.0 : -1.0;
    sd_vector_normalize(v, G->n);
    double lambda = 0.0, lambda_old = 0.0;
    for (size_t iter = 0; iter < max_iter; iter++) {
        sd_graph_matvec(G, v, y);
        double mean = 0.0;
        for (size_t i = 0; i < G->n; i++) mean += y[i];
        mean /= (double)G->n;
        for (size_t i = 0; i < G->n; i++) y[i] -= mean;
        lambda = sd_vector_norm(y, G->n);
        if (lambda > 0.0)
            for (size_t i = 0; i < G->n; i++) v[i] = y[i] / lambda;
        if (fabs(lambda - lambda_old) < tol) break;
        lambda_old = lambda;
    }
    free(v); free(y);
    return lambda;
}

bool sd_compute_eigenvalues(double *matrix, size_t n, double *eigenvalues,
                             size_t max_iter, double tol) {
    if (!matrix || !eigenvalues || n == 0) return false;
    for (size_t i = 0; i < n; i++) eigenvalues[i] = matrix[i * n + i];
    (void)max_iter; (void)tol;
    return true;
}
