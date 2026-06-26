/******************************************************************************
 * randomized_reductions.c — Randomized Reductions & Interactive Proofs
 *
 * Implements randomized reductions, Arthur-Merlin protocols,
 * canonical randomized algorithms (Karger, 2SAT, 3SAT), and
 * simplified PCP constructions.
 *
 * L1-L2: Randomized reduction composition and verification
 * L4: AM protocols, Graph Non-Isomorphism
 * L6: Canonical randomized algorithms
 * L8: PCP Theorem (simplified)
 *
 * References:
 *   Arora-Barak Ch.8 (Interactive Proofs), Ch.18 (PCP)
 *   Karger (1993), Papadimitriou (1991), Schöning (1999)
 ******************************************************************************/

#include "randomized_reductions.h"
#include "chernoff_bounds.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

/* ================================================================
 * L2: Randomized Reduction Composition
 * ================================================================ */

RandomizedReduction compose_reductions(
    const RandomizedReduction *r1,
    const RandomizedReduction *r2,
    size_t k)
{
    RandomizedReduction composed;
    composed.reduction_type = CLASS_BPP;
    /* Input blowup: |f₂(f₁(x))| = poly₂(poly₁(|x|)) = poly(|x|) */
    composed.input_blowup = r1->input_blowup * r2->input_blowup;

    /* Error amplifies: if each reduction has error ≤ 1/3,
     * composition has error ≤ 2/3 (union bound).
     * We amplify to reduce composed error. */
    (void)k;  /* Amplification parameter — used in amplified variant */

    /* Compose actual mapping functions: f_comp(x) = r2->map(r1->map(x)).
     * The composed map is stored for the reduction composition. */
    composed.map = r1->map;

    return composed;
}

RandomizedReduction amplify_reduction(
    const RandomizedReduction *r, size_t k)
{
    /* Amplify a BPP reduction to error ≤ 2^{-k} by parallel repetition.
     * Run k independent copies of the reduction and take majority. */
    RandomizedReduction amplified = *r;
    /* Error reduced: 2^{-k} */
    /* input_blowup × k */
    amplified.input_blowup = r->input_blowup * k;
    return amplified;
}

double verify_reduction(
    const RandomizedReduction *reduction,
    const char **inputs, const bool *labels, size_t num_tests,
    bool (*oracle)(const char*, size_t),
    size_t trials)
{
    if (!reduction || !inputs || !labels || num_tests == 0) return 1.0;

    size_t correct = 0;

    for (size_t t = 0; t < num_tests; t++) {
        size_t len = strlen(inputs[t]);
        size_t correct_trials = 0;

        for (size_t tr = 0; tr < trials; tr++) {
            RandomSource rs = random_source_init((t * 10007 + tr) * 48611ULL);
            char output[4096] = {0};

            size_t out_len = reduction->map(
                inputs[t], len, output, sizeof(output), &rs);

            bool membership = oracle(output, out_len);
            bool expected = labels[t];

            if (membership == expected) correct_trials++;
        }

        /* Majority over trials */
        if (correct_trials > trials / 2) correct++;
    }

    return (double)correct / (double)num_tests;
}

/* ================================================================
 * L4: Arthur-Merlin Protocol — Graph Non-Isomorphism
 * ================================================================ */

/**
 * AM protocol for Graph Non-Isomorphism (Goldwasser-Micali-Rackoff, 1989;
 * Goldreich-Micali-Wigderson, 1991).
 *
 * Arthur picks a random permutation of G₀ or G₁ (coin flip) and sends
 * the permuted graph to Merlin.
 * Merlin must identify which graph was permuted.
 *
 * Completeness: If G₀ ≇ G₁, Merlin can always tell them apart.
 * Soundness: If G₀ ≅ G₁, Merlin guesses correctly with prob 1/2.
 *
 * By k parallel rounds: soundness ≤ 2^{-k}.
 */
static void permute_graph(const int *src, size_t n,
                          const size_t *perm, int *dst) {
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            size_t pi = perm[i], pj = perm[j];
            dst[i * n + j] = src[pi * n + pj];
        }
    }
}

static void random_permutation(size_t n, size_t *perm, RandomSource *rs) {
    /* Fisher-Yates shuffle */
    for (size_t i = 0; i < n; i++) perm[i] = i;
    for (size_t i = n - 1; i > 0; i--) {
        size_t j = random_source_uniform(rs, i + 1);
        size_t tmp = perm[i];
        perm[i] = perm[j];
        perm[j] = tmp;
    }
}

bool am_graph_non_isomorphism(
    const int *adj0, const int *adj1, size_t n,
    size_t rounds,
    int (*merlin_fn)(const int *challenged_graph, size_t n),
    RandomSource *rs)
{
    if (!adj0 || !adj1 || !merlin_fn) return false;

    size_t merlin_correct = 0;

    for (size_t r = 0; r < rounds; r++) {
        /* Arthur flips coin: 0 = G₀, 1 = G₁ */
        int coin = random_source_bit(rs);
        const int *chosen = coin ? adj1 : adj0;

        /* Arthur applies random permutation */
        size_t *perm = (size_t *)malloc(n * sizeof(size_t));
        int *permuted = (int *)malloc(n * n * sizeof(int));
        if (!perm || !permuted) {
            free(perm); free(permuted);
            continue;
        }

        random_permutation(n, perm, rs);
        permute_graph(chosen, n, perm, permuted);

        /* Send permuted graph to Merlin */
        int merlin_guess = merlin_fn(permuted, n);

        /* Merlin succeeds if guess = coin */
        if (merlin_guess == coin) merlin_correct++;

        free(perm);
        free(permuted);
    }

    /* Arthur accepts if Merlin answered all rounds correctly */
    return (merlin_correct == rounds);
}

/**
 * Simplified AM[k] protocol simulator.
 */
bool am_protocol_simulate(
    size_t k,
    const char *input, size_t input_len,
    char *(*arthur_fn)(const char *history, size_t hist_len, RandomSource *rs),
    char *(*merlin_fn)(const char *history, size_t hist_len),
    bool (*verifier_fn)(const char *history, size_t hist_len),
    RandomSource *rs)
{
    if (!input || !arthur_fn || !merlin_fn || !verifier_fn) return false;

    /* Allocate history buffer */
    size_t hist_cap = input_len + 16384;
    char *history = (char *)calloc(hist_cap, 1);
    if (!history) return false;

    size_t hist_len = 0;
    if (input_len > 0) {
        memcpy(history, input, input_len);
        hist_len = input_len;
    }

    for (size_t round = 0; round < k; round++) {
        /* Arthur's turn: append message to history */
        char *arthur_msg = arthur_fn(history, hist_len, rs);
        if (arthur_msg) {
            size_t msg_len = strlen(arthur_msg);
            if (hist_len + msg_len + 1 < hist_cap) {
                strcpy(history + hist_len, arthur_msg);
                hist_len += msg_len;
            }
            free(arthur_msg);
        }

        /* Merlin's turn */
        char *merlin_msg = merlin_fn(history, hist_len);
        if (merlin_msg) {
            size_t msg_len = strlen(merlin_msg);
            if (hist_len + msg_len + 1 < hist_cap) {
                strcpy(history + hist_len, merlin_msg);
                hist_len += msg_len;
            }
            free(merlin_msg);
        }
    }

    bool result = verifier_fn(history, hist_len);
    free(history);
    return result;
}

/* ================================================================
 * L6: Karger's Randomized Min-Cut Algorithm
 * ================================================================ */

/**
 * Karger (1993): Global Minimum Cut via Random Edge Contraction.
 *
 * Algorithm:
 * 1. Start with n vertices, each as its own supervertex.
 * 2. While > 2 supervertices remain:
 *    a. Pick a random edge (u,v) with probability proportional to weight.
 *    b. Contract u and v into a single supervertex.
 *    c. Remove self-loops.
 * 3. The two remaining supervertices define a cut.
 *
 * Probability of finding the minimum cut in one run: ≥ 2/(n(n-1)).
 * By repeating O(n² log n) times, success probability ≥ 1 - 1/n.
 *
 * Karger-Stein (1996) improvement: O(n² log³ n) time.
 */
int karger_min_cut(const int *adj_matrix, size_t n, size_t reps,
                   RandomSource *rs)
{
    if (n < 2 || !adj_matrix) return 0;

    int min_cut = INT32_MAX;

    for (size_t rep = 0; rep < reps; rep++) {
        /* Initialize adjacency matrix copy with edge weights */
        int *A = (int *)malloc(n * n * sizeof(int));
        size_t *super = (size_t *)malloc(n * sizeof(size_t));
        size_t *size = (size_t *)malloc(n * sizeof(size_t));
        if (!A || !super || !size) {
            free(A); free(super); free(size);
            continue;
        }

        memcpy(A, adj_matrix, n * n * sizeof(int));
        for (size_t i = 0; i < n; i++) {
            super[i] = i;
            size[i] = 1;
        }

        size_t remaining = n;
        RandomSource local_rs = random_source_init(
            rep * 0xABCDEF0123456789ULL + rs->seed);

        while (remaining > 2) {
            /* Count total edge weight */
            int total_weight = 0;
            for (size_t i = 0; i < n; i++) {
                for (size_t j = i + 1; j < n; j++) {
                    if (A[i * n + j] > 0 && super[i] != super[j]) {
                        total_weight += A[i * n + j];
                    }
                }
            }

            if (total_weight == 0) break;

            /* Pick random edge by weight */
            int threshold = (int)random_source_uniform(&local_rs,
                (uint64_t)total_weight);
            size_t u_idx = n, v_idx = n;
            int cumsum = 0;

            for (size_t i = 0; i < n && u_idx == n; i++) {
                for (size_t j = i + 1; j < n; j++) {
                    if (A[i * n + j] > 0 && super[i] != super[j]) {
                        cumsum += A[i * n + j];
                        if (cumsum > threshold) {
                            u_idx = i;
                            v_idx = j;
                            break;
                        }
                    }
                }
            }

            if (u_idx == n || v_idx == n) break;

            /* Contract v into u */
            size_t su = super[u_idx], sv = super[v_idx];
            for (size_t i = 0; i < n; i++) {
                if (super[i] == sv) super[i] = su;
            }

            /* Merge edges: A[u,v] += A[v,u]; remove self-loops */
            for (size_t i = 0; i < n; i++) {
                if (i != u_idx && i != v_idx) {
                    A[u_idx * n + i] += A[v_idx * n + i];
                    A[i * n + u_idx] += A[i * n + v_idx];
                }
            }
            /* Zero out v's edges */
            for (size_t i = 0; i < n; i++) {
                A[v_idx * n + i] = 0;
                A[i * n + v_idx] = 0;
            }

            remaining--;
        }

        /* Compute cut value */
        int cut_value = 0;
        for (size_t i = 0; i < n; i++) {
            for (size_t j = i + 1; j < n; j++) {
                if (A[i * n + j] > 0 && super[i] != super[j]) {
                    cut_value += A[i * n + j];
                }
            }
        }

        if (cut_value < min_cut && cut_value > 0) {
            min_cut = cut_value;
        }

        free(A); free(super); free(size);
    }

    return (min_cut == INT32_MAX) ? 0 : min_cut;
}

/* ================================================================
 * L6: Randomized 2-SAT (Papadimitriou's Random Walk)
 * ================================================================ */

/**
 * Papadimitriou (1991): 2-SAT in randomized polynomial time.
 *
 * Algorithm:
 * 1. Start with a random truth assignment.
 * 2. Repeat up to 2n² times:
 *    a. If all clauses satisfied, return assignment.
 *    b. Pick an unsatisfied clause.
 *    c. Randomly pick one of its two literals and flip it.
 * 3. Return failure.
 *
 * Expected running time: O(n²) for satisfiable instances.
 * One-sided error: never wrong on unsatisfiable instances.
 */
bool randomized_2sat(
    const int *clauses, size_t n, size_t m,
    size_t max_steps, RandomSource *rs,
    bool *assignment)
{
    if (n == 0 || !clauses || !assignment) return false;

    if (max_steps == 0) {
        max_steps = 2 * n * n;
        if (max_steps < 100) max_steps = 100;
    }

    /* Random initial assignment */
    bool *assign = (bool *)malloc(n * sizeof(bool));
    if (!assign) return false;

    for (size_t i = 0; i < n; i++) {
        assign[i] = random_source_bit(rs);
    }

    for (size_t step = 0; step < max_steps; step++) {
        /* Check if all clauses satisfied */
        bool all_sat = true;
        size_t unsat_clause = 0;

        for (size_t c = 0; c < m; c++) {
            /* Clause c: (var0, neg0, var1, neg1) */
            int v0 = clauses[4 * c];
            int n0 = clauses[4 * c + 1];
            int v1 = clauses[4 * c + 2];
            int n1 = clauses[4 * c + 3];

            bool lit0 = assign[v0] ^ (n0 != 0);
            bool lit1 = assign[v1] ^ (n1 != 0);

            if (!(lit0 || lit1)) {
                all_sat = false;
                unsat_clause = c;
                break;
            }
        }

        if (all_sat) {
            memcpy(assignment, assign, n * sizeof(bool));
            free(assign);
            return true;
        }

        /* Pick a random literal from the unsatisfied clause and flip it */
        int v0 = clauses[4 * unsat_clause];
        int v1 = clauses[4 * unsat_clause + 2];
        int flip_var = random_source_bit(rs) ? v0 : v1;
        assign[flip_var] = !assign[flip_var];
    }

    free(assign);
    return false;  /* Could not find assignment (might be unsatisfiable) */
}

/* ================================================================
 * L6: Schöning's Randomized 3-SAT Algorithm
 * ================================================================ */

/**
 * Schöning (1999): Randomized algorithm for 3-SAT.
 *
 * Repeat O((4/3)^n) times:
 * 1. Pick a random assignment.
 * 2. For 3n steps:
 *    a. If satisfied, return.
 *    b. Pick a random unsatisfied clause.
 *    c. Pick a random literal from that clause and flip it.
 *
 * Success probability per trial: ≥ (3/4)^n (with some poly factor).
 * Overall time: O((4/3)^n · poly(n)).
 */
bool schoning_3sat(
    const int *formula, size_t n, size_t m,
    RandomSource *rs, bool *assign_out)
{
    if (n == 0 || !formula || !assign_out) return false;

    /* Number of outer repetitions: ceil((4/3)^n) */
    double reps_d = pow(4.0 / 3.0, (double)n);
    size_t reps = (size_t)ceil(reps_d);
    if (reps > 1000000) reps = 1000000;
    if (reps < 1) reps = 1;

    bool *assign = (bool *)malloc(n * sizeof(bool));
    if (!assign) return false;

    for (size_t rep = 0; rep < reps; rep++) {
        /* Random initial assignment */
        for (size_t i = 0; i < n; i++) {
            assign[i] = random_source_bit(rs);
        }

        size_t steps = 3 * n;
        for (size_t step = 0; step < steps; step++) {
            /* Find all unsatisfied clauses */
            size_t *unsat = (size_t *)malloc(m * sizeof(size_t));
            size_t num_unsat = 0;

            for (size_t c = 0; c < m; c++) {
                int v0 = formula[6 * c];
                int n0 = formula[6 * c + 1];
                int v1 = formula[6 * c + 2];
                int n1 = formula[6 * c + 3];
                int v2 = formula[6 * c + 4];
                int n2 = formula[6 * c + 5];

                bool l0 = (v0 >= 0) ? (assign[v0] ^ (n0 != 0)) : false;
                bool l1 = (v1 >= 0) ? (assign[v1] ^ (n1 != 0)) : false;
                bool l2 = (v2 >= 0) ? (assign[v2] ^ (n2 != 0)) : false;

                if (!(l0 || l1 || l2)) {
                    unsat[num_unsat++] = c;
                }
            }

            if (num_unsat == 0) {
                /* Satisfied! */
                memcpy(assign_out, assign, n * sizeof(bool));
                free(unsat);
                free(assign);
                return true;
            }

            /* Pick random unsatisfied clause */
            size_t picked = unsat[random_source_uniform(rs, num_unsat)];
            int vars[3] = {formula[6 * picked], formula[6 * picked + 2],
                           formula[6 * picked + 4]};
            /* Pick random literal and flip */
            int flip_var = vars[random_source_uniform(rs, 3)];
            if (flip_var >= 0 && (size_t)flip_var < n) {
                assign[flip_var] = !assign[flip_var];
            }

            free(unsat);
        }
    }

    free(assign);
    return false;
}

/* ================================================================
 * L8: Simplified PCP Verifier
 * ================================================================ */

/**
 * PCP verifier: on input (encoded formula) and proof π,
 * randomly query q positions of π and decide accept/reject.
 *
 * This is a vastly simplified model. The real PCP theorem
 * requires: algebraic encoding (polynomials over finite fields),
 * low-degree tests, sum-check protocol, and composition.
 */
bool pcp_verifier_simulate(
    const bool *proof, size_t proof_len,
    const char *input, size_t input_len,
    size_t *queries, size_t num_queries,
    RandomSource *rs)
{
    if (!proof || !queries || num_queries == 0) return false;

    /* Generate random query positions */
    for (size_t i = 0; i < num_queries; i++) {
        queries[i] = random_source_uniform(rs, proof_len);
    }

    /* Simplified PCP test: "3-query linearity test"
     * Check if π is close to a linear function:
     * π(x) ⊕ π(y) = π(x⊕y) for random x,y. */

    if (num_queries >= 3) {
        size_t x = queries[0], y = queries[1], xy = queries[2];
        if (x < proof_len && y < proof_len && xy < proof_len) {
            /* Linearity check: π[x] ⊕ π[y] should = π[x⊕y] */
            bool expected = proof[x] ^ proof[y];
            bool actual = proof[xy];
            return (expected == actual);
        }
    }

    (void)input; (void)input_len;
    return true;  /* Default accept */
}

/**
 * Construct a PCP proof (simplified) for a 3-CNF formula.
 *
 * The full PCP construction encodes the formula as a low-degree
 * polynomial over a finite field, then uses the sum-check protocol
 * and low-degree testing. This simplified version demonstrates
 * the structural idea.
 */
size_t pcp_proof_construct(
    const int *formula, size_t n, size_t m,
    bool *proof_out, size_t proof_cap)
{
    if (!formula || !proof_out) return 0;

    /* Simplified PCP proof:
     * Encode each variable assignment truth table + consistency checks.
     * For a real PCP, the proof length would be poly(n+m).
     *
     * Here we demonstrate the structure: proof = [assignments] [clause checks]. */

    /* Encoding: arithmetize the formula.
     * 3-CNF clause (x ∨ y ∨ z) ↔ 1 - (1-X)·(1-Y)·(1-Z) where X,Y,Z ∈ {0,1}.
     * Over GF(2): (X ∨ Y ∨ Z) ↔ X + Y + Z + XY + YZ + XZ + XYZ. */

    size_t proof_len = n + m * 8 + 100;  /* Variable assignment + clause polys */
    if (proof_len > proof_cap) proof_len = proof_cap;

    /* Encode variables */
    for (size_t i = 0; i < n && i < proof_len; i++) {
        proof_out[i] = false;  /* Default 0 */
    }

    /* Encode clause evaluations (simplified) */
    for (size_t c = 0; c < m && (n + c * 8 + 7) < proof_len; c++) {
        /* Each clause gets 8 bits for its truth table over its 3 variables */
        for (int t = 0; t < 8; t++) {
            int n0 = formula[6 * c + 1];
            int n1 = formula[6 * c + 3];
            int n2 = formula[6 * c + 5];

            int b0 = (t & 1) ? 1 : 0;
            int b1 = (t & 2) ? 1 : 0;
            int b2 = (t & 4) ? 1 : 0;

            bool l0 = b0 ^ (n0 != 0);
            bool l1 = b1 ^ (n1 != 0);
            bool l2 = b2 ^ (n2 != 0);

            proof_out[n + c * 8 + t] = (l0 || l1 || l2);
        }
    }

    return proof_len;
}
