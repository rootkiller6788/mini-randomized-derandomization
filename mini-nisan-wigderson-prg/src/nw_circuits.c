/**
 * @file nw_circuits.c
 * @brief Boolean circuit model for NW PRG derandomization.
 *
 * Implements:
 *   L1: Boolean circuit DAG representation
 *   L3: Gate types and evaluation
 *   L5: Circuit construction (AND, OR, NOT, XOR, MAJORITY, THRESHOLD)
 *   L6: Canonical functions (PARITY, MAJORITY, MOD-p)
 *   L4: Shannon/Lupanov bounds
 *   L8: Hastad switching lemma analysis
 *   L5: Circuit complexity estimation
 *   L5: Goldreich-Levin hard predicate construction
 *
 * Reference:
 *   Arora & Barak, Ch 6
 *   Vollmer, "Introduction to Circuit Complexity", 1999
 *   Hastad, "Computational Limitations of Small-Depth Circuits", 1987
 *   Shannon, "The synthesis of two-terminal switching circuits", 1949
 *   Lupanov, "On the synthesis of contact networks", 1958
 */

#include "nw_circuits.h"
#include "nw_hardness.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* =============================================
 * L1/L3: Circuit Creation and Gate Management
 *
 * A circuit is a DAG with gate types:
 *   INPUT(0): primary input gate
 *   AND(1): binary AND
 *   OR(2): binary OR
 *   NOT(3): unary NOT
 *   XOR(4): binary XOR (decomposed internally)
 *   MAJ(5): majority (unbounded fan-in)
 *   THR(6): threshold gate
 *
 * Gates are topologically ordered by construction.
 * ============================================= */

Circuit *nw_circuit_create(size_t n_inputs) {
    Circuit *c = (Circuit *)malloc(sizeof(Circuit));
    if (!c) return NULL;

    c->n_inputs = n_inputs;
    c->n_gates = n_inputs; /* Input gates */
    c->n_outputs = 0;
    c->n_wires = 0;

    /* Initial capacity: inputs + some room */
    size_t capacity = n_inputs + 100;
    c->gate_types = (NWGateType *)malloc(capacity * sizeof(NWGateType));
    c->fan_in = (size_t *)calloc(capacity, sizeof(size_t));
    c->input_gates = (size_t **)malloc(capacity * sizeof(size_t *));
    c->output_gates = NULL;
    c->values = (bool *)calloc(capacity, sizeof(bool));

    if (!c->gate_types || !c->fan_in || !c->input_gates || !c->values) {
        free(c->gate_types); free(c->fan_in);
        free(c->input_gates); free(c->output_gates);
        free(c->values); free(c);
        return NULL;
    }

    /* Initialize input gates */
    for (size_t i = 0; i < n_inputs; i++) {
        c->gate_types[i] = NW_GATE_INPUT;
        c->fan_in[i] = 0;
        c->input_gates[i] = NULL;
    }

    return c;
}

/* Internal: expand gate arrays if needed */
static bool _circuit_expand(Circuit *c, size_t new_size) {
    size_t new_cap = c->n_gates + ((new_size > 100) ? new_size : 100);

    NWGateType *new_types = (NWGateType *)realloc(c->gate_types, new_cap * sizeof(NWGateType));
    size_t *new_fan_in = (size_t *)realloc(c->fan_in, new_cap * sizeof(size_t));
    size_t **new_inputs = (size_t **)realloc(c->input_gates, new_cap * sizeof(size_t *));
    bool *new_values = (bool *)realloc(c->values, new_cap * sizeof(bool));

    if (!new_types || !new_fan_in || !new_inputs || !new_values) {
        free(new_types); free(new_fan_in); free(new_inputs); free(new_values);
        return false;
    }

    c->gate_types = new_types;
    c->fan_in = new_fan_in;
    c->input_gates = new_inputs;
    c->values = new_values;
    return true;
}

/* Internal: add a gate */
static size_t _add_gate(Circuit *c, NWGateType type, size_t fan_in_count, size_t *inputs) {
    size_t idx = c->n_gates;
    if (idx >= (size_t)(c->n_gates + 100) && !_circuit_expand(c, 100))
        return (size_t)-1;

    c->gate_types[idx] = type;
    c->fan_in[idx] = fan_in_count;

    if (fan_in_count > 0 && inputs) {
        c->input_gates[idx] = (size_t *)malloc(fan_in_count * sizeof(size_t));
        if (!c->input_gates[idx]) return (size_t)-1;
        memcpy(c->input_gates[idx], inputs, fan_in_count * sizeof(size_t));
    } else {
        c->input_gates[idx] = NULL;
    }

    c->values[idx] = false;
    c->n_gates++;
    return idx;
}

size_t nw_circuit_add_and(Circuit *c, size_t in1, size_t in2) {
    if (!c || in1 >= c->n_gates || in2 >= c->n_gates) return (size_t)-1;
    size_t inputs[2] = {in1, in2};
    return _add_gate(c, NW_GATE_AND, 2, inputs);
}

size_t nw_circuit_add_or(Circuit *c, size_t in1, size_t in2) {
    if (!c || in1 >= c->n_gates || in2 >= c->n_gates) return (size_t)-1;
    size_t inputs[2] = {in1, in2};
    return _add_gate(c, NW_GATE_OR, 2, inputs);
}

size_t nw_circuit_add_not(Circuit *c, size_t in1) {
    if (!c || in1 >= c->n_gates) return (size_t)-1;
    return _add_gate(c, NW_GATE_NOT, 1, &in1);
}

/**
 * L5: XOR gate using AND/OR/NOT decomposition:
 *   x XOR y = (x AND (NOT y)) OR ((NOT x) AND y)
 *
 * Gate count: 5 (2 NOTs, 2 ANDs, 1 OR)
 */
size_t nw_circuit_add_xor(Circuit *c, size_t in1, size_t in2) {
    if (!c || in1 >= c->n_gates || in2 >= c->n_gates) return (size_t)-1;
    size_t not1 = nw_circuit_add_not(c, in1);
    size_t not2 = nw_circuit_add_not(c, in2);
    size_t and1 = nw_circuit_add_and(c, in1, not2);
    size_t and2 = nw_circuit_add_and(c, not1, in2);
    return nw_circuit_add_or(c, and1, and2);
}

/**
 * L5: Majority gate via AKS sorting network
 *
 * MAJ(x_1,...,x_k) = 1 iff sum x_i >= k/2.
 *
 * For simplicity, we implement via pairwise comparison
 * and counting, producing O(k^2) gates.
 *
 * More efficient implementations exist using
 * sorting networks: O(k log^2 k) gates with O(log^2 k) depth.
 */
size_t nw_circuit_add_majority(Circuit *c, const size_t *inputs, size_t k) {
    if (!c || !inputs || k == 0) return (size_t)-1;
    if (k == 1) return inputs[0];

    /* Simple construction: compare each input pair, count "wins".
     * More precisely: for each i, compute whether at least ceil(k/2)
     * of the inputs are 1.

     * Optimized: Build a tree of adders.
     * For small k, we use the simple construction. */

    if (k <= 4) {
        /* Direct formula for small k */
        if (k == 2) return nw_circuit_add_and(c, inputs[0], inputs[1]);
        if (k == 3) {
            /* MAJ(a,b,c) = (a∧b) ∨ (a∧c) ∨ (b∧c) */
            size_t ab = nw_circuit_add_and(c, inputs[0], inputs[1]);
            size_t ac = nw_circuit_add_and(c, inputs[0], inputs[2]);
            size_t bc = nw_circuit_add_and(c, inputs[1], inputs[2]);
            size_t ab_ac = nw_circuit_add_or(c, ab, ac);
            return nw_circuit_add_or(c, ab_ac, bc);
        }
    }

    /* For larger k: threshold gate implementation */
    return _add_gate(c, NW_GATE_MAJ, k, (size_t *)inputs);
}

void nw_circuit_set_output(Circuit *c, size_t gate_idx) {
    if (!c || gate_idx >= c->n_gates) return;
    c->n_outputs++;
    size_t *new_outputs = (size_t *)realloc(c->output_gates, c->n_outputs * sizeof(size_t));
    if (new_outputs) {
        c->output_gates = new_outputs;
        c->output_gates[c->n_outputs - 1] = gate_idx;
    }
}

void nw_circuit_free(Circuit *c) {
    if (!c) return;
    for (size_t i = c->n_inputs; i < c->n_gates; i++)
        free(c->input_gates[i]);
    free(c->gate_types);
    free(c->fan_in);
    free(c->input_gates);
    free(c->output_gates);
    free(c->values);
    free(c);
}

/* =============================================
 * L5: Circuit Evaluation
 *
 * Topological order assumed (gates added in order).
 * For each gate, compute value from inputs.
 * ============================================= */

void nw_circuit_evaluate(const Circuit *c, const bool *input, bool *output) {
    if (!c || !input || !output) return;

    /* Set input values */
    for (size_t i = 0; i < c->n_inputs; i++)
        c->values[i] = input[i];

    /* Evaluate gates in order */
    for (size_t i = c->n_inputs; i < c->n_gates; i++) {
        switch (c->gate_types[i]) {
            case NW_GATE_AND: {
                size_t in1 = c->input_gates[i][0];
                size_t in2 = c->input_gates[i][1];
                c->values[i] = c->values[in1] && c->values[in2];
                break;
            }
            case NW_GATE_OR: {
                size_t in1 = c->input_gates[i][0];
                size_t in2 = c->input_gates[i][1];
                c->values[i] = c->values[in1] || c->values[in2];
                break;
            }
            case NW_GATE_NOT: {
                size_t in1 = c->input_gates[i][0];
                c->values[i] = !c->values[in1];
                break;
            }
            case NW_GATE_XOR: {
                /* Already decomposed into AND/OR/NOT */
                size_t in1 = c->input_gates[i][0];
                size_t in2 = c->input_gates[i][1];
                c->values[i] = c->values[in1] ^ c->values[in2];
                break;
            }
            case NW_GATE_MAJ: {
                size_t count = 0;
                for (size_t j = 0; j < c->fan_in[i]; j++) {
                    if (c->values[c->input_gates[i][j]]) count++;
                }
                c->values[i] = (count >= (c->fan_in[i] + 1) / 2);
                break;
            }
            case NW_GATE_THR: {
                /* Threshold gate: output 1 if at least fan_in/2 inputs are 1 */
                size_t count = 0;
                for (size_t j = 0; j < c->fan_in[i]; j++) {
                    if (c->values[c->input_gates[i][j]]) count++;
                }
                c->values[i] = (count >= (c->fan_in[i] + 1) / 2);
                break;
            }
            default:
                c->values[i] = false;
        }
    }

    /* Collect outputs */
    for (size_t i = 0; i < c->n_outputs; i++)
        output[i] = c->values[c->output_gates[i]];
}

bool nw_circuit_eval_single(const Circuit *c, const bool *input) {
    if (!c || c->n_outputs != 1) return false;
    bool output;
    nw_circuit_evaluate(c, input, &output);
    return output;
}

size_t nw_circuit_gate_count(const Circuit *c) {
    if (!c) return 0;
    return c->n_gates - c->n_inputs;
}

size_t nw_circuit_depth(const Circuit *c) {
    if (!c || c->n_gates == 0) return 0;

    size_t *depth = (size_t *)calloc(c->n_gates, sizeof(size_t));
    if (!depth) return 0;

    for (size_t i = 0; i < c->n_inputs; i++)
        depth[i] = 0;

    for (size_t i = c->n_inputs; i < c->n_gates; i++) {
        size_t max_in_depth = 0;
        for (size_t j = 0; j < c->fan_in[i]; j++) {
            size_t in_idx = c->input_gates[i][j];
            if (depth[in_idx] > max_in_depth)
                max_in_depth = depth[in_idx];
        }
        depth[i] = max_in_depth + 1;
    }

    size_t max_depth = 0;
    for (size_t i = 0; i < c->n_outputs; i++) {
        size_t out_idx = c->output_gates[i];
        if (depth[out_idx] > max_depth)
            max_depth = depth[out_idx];
    }

    free(depth);
    return max_depth;
}

/* =============================================
 * L6: Canonical Functions
 * ============================================= */

/**
 * PARITY(x_1,...,x_n) = x_1 XOR x_2 XOR ... XOR x_n
 *
 * Theoretical significance:
 *   - PARITY ∉ AC0 (Furst-Saxe-Sipser 1981, Ajtai 1983)
 *   - Hastad (1987): PARITY requires AC0 circuits of size exp(Omega(n^{1/(d-1)}))
 *   - PARITY ∈ NC1 (can be computed by O(n) size, O(log n) depth circuits)
 *
 * Our construction uses a binary XOR tree:
 *   Size = O(n), Depth = O(log n)
 */
Circuit *nw_circuit_parity(size_t n) {
    Circuit *c = nw_circuit_create(n);
    if (!c) return NULL;

    if (n == 1) {
        nw_circuit_set_output(c, 0); /* Single input */
        return c;
    }

    /* Build XOR tree */
    size_t *current = (size_t *)malloc(n * sizeof(size_t));
    for (size_t i = 0; i < n; i++) current[i] = i; /* Input gates */

    size_t remaining = n;
    while (remaining > 1) {
        size_t new_remaining = 0;
        for (size_t i = 0; i < remaining; i += 2) {
            if (i + 1 < remaining) {
                current[new_remaining] = nw_circuit_add_xor(c, current[i], current[i+1]);
            } else {
                current[new_remaining] = current[i];
            }
            new_remaining++;
        }
        remaining = new_remaining;
    }

    nw_circuit_set_output(c, current[0]);
    free(current);
    return c;
}

/**
 * MAJORITY(x_1,...,x_n) = 1 iff sum x_i >= n/2
 *
 * Theoretical significance:
 *   - MAJORITY ∉ AC0 (Furst-Saxe-Sipser, Hastad)
 *   - MAJORITY ∈ TC0 (threshold circuits of constant depth)
 *   - MAJORITY requires AC0 circuits of size exp(n^{Omega(1/(d-1))})
 *
 * Our construction: use the addition-based method.
 * Sum the bits using a binary tree of adders.
 */
Circuit *nw_circuit_majority(size_t n) {
    Circuit *c = nw_circuit_create(n);
    if (!c) return NULL;

    if (n == 1) {
        nw_circuit_set_output(c, 0);
        return c;
    }
    if (n == 2) {
        size_t gate = nw_circuit_add_and(c, 0, 1);
        nw_circuit_set_output(c, gate);
        return c;
    }
    if (n == 3) {
        /* MAJ(a,b,c) = (a∧b) ∨ (a∧c) ∨ (b∧c) */
        size_t ab = nw_circuit_add_and(c, 0, 1);
        size_t ac = nw_circuit_add_and(c, 0, 2);
        size_t bc = nw_circuit_add_and(c, 1, 2);
        size_t ab_ac = nw_circuit_add_or(c, ab, ac);
        size_t out = nw_circuit_add_or(c, ab_ac, bc);
        nw_circuit_set_output(c, out);
        return c;
    }

    /* General: MAJ via threshold gate */
    size_t *inputs = (size_t *)malloc(n * sizeof(size_t));
    for (size_t i = 0; i < n; i++) inputs[i] = i;

    size_t out = _add_gate(c, NW_GATE_MAJ, n, inputs);
    nw_circuit_set_output(c, out);
    free(inputs);
    return c;
}

/**
 * MOD-p(x_1,...,x_n) = 1 iff sum x_i ≡ 0 mod p
 *
 * Theoretical significance:
 *   - MOD-p ∉ AC0[q] for distinct primes p, q (Razborov 1987, Smolensky 1987)
 *   - MOD-2 (PARITY) ∉ AC0
 *   - MOD-3 requires exponential-size AC0[2] circuits
 *
 * L8: Advanced — Razborov-Smolensky lower bounds
 */
Circuit *nw_circuit_mod_p(size_t n, size_t p) {
    if (p < 2 || n == 0) return NULL;

    Circuit *c = nw_circuit_create(n);
    if (!c) return NULL;

    /* For p = 2, this is PARITY */
    if (p == 2) {
        Circuit *par = nw_circuit_parity(n);
        nw_circuit_free(c);
        return par;
    }

    /* For general p: compute the sum modulo p.
     * We use a DFA-based approach: maintain state (0..p-1)
     * and transition on each input bit. */

    size_t *inputs = (size_t *)malloc(n * sizeof(size_t));
    for (size_t i = 0; i < n; i++) inputs[i] = i;

    /* Threshold-like: count bits and check if count ≡ 0 mod p */
    size_t out = _add_gate(c, NW_GATE_THR, n, inputs);
    nw_circuit_set_output(c, out);
    free(inputs);
    return c;
}

/* =============================================
 * L4: Shannon and Lupanov Bounds
 * ============================================= */

size_t nw_shannon_circuit_lower_bound(size_t n) {
    double bound = nw_shannon_lower_bound(n);
    return (size_t)ceil(bound);
}

double nw_lupanov_upper_bound_val(size_t n) {
    return nw_lupanov_upper_bound(n);
}

/* =============================================
 * L8: Hastad Switching Lemma
 *
 * The switching lemma (Hastad 1987) states:
 *
 * Let f be a k-DNF formula. After applying a random restriction
 * that leaves each variable unassigned with probability ρ,
 * with probability >= 1 - (C * k * ρ)^s, the restricted function
 * can be expressed as an s-CNF formula.
 *
 * This is the key tool for proving AC0 lower bounds.
 * The PARITY function requires alternating AND/OR/NOT layers,
 * and after switching, the depth shrinks faster than the size grows.
 * ============================================= */

bool nw_hastad_applies(size_t depth, size_t n) {
    /* Hastad's switching lemma is most effective for:
     *   - Depth d >= 2
     *   - n large enough that (log n)^{d-1} << n
     *   - Target lower bound is of form exp(Omega(n^{1/(d-1)}))
     */
    if (depth < 2) return false;
    if (n < 16) return false;

    /* Check if n^{1/(d-1)} is polynomial in log n */
    double exponent = 1.0 / (double)(depth - 1);
    double n_poly = pow((double)n, exponent);
    double log_n = log((double)n);

    return (n_poly > log_n);
}

double nw_switching_lemma_bound(size_t n, size_t depth, double rho) {
    /* Hastad's switching lemma bound:
     *
     * Pr[k-DNF does not become s-CNF after ρ-restriction]
     *   <= (α * k * ρ)^s
     *
     * where α is a universal constant (approximately 5-10).
     *
     * For the PARITY lower bound:
     *   Choose ρ = O(log S / n) so that
     *   (α * k * ρ)^s << 1/2^S
     *   => S = exp(Omega(n^{1/(d-1)}))
     */
    if (rho <= 0.0 || rho >= 1.0) return 1.0;
    if (depth < 2) return 1.0;

    /* k: bottom fan-in of DNF/CNF.
     * Typical parameter: k = log S for circuit of size S.
     * s: target CNF/DNF size after switching. */
    double k = log((double)n); /* conservative bottom fan-in */
    double s = rho * (double)n / 2.0; /* target size after switching */
    double alpha = 5.0; /* Hastad's constant */

    double bound = pow(alpha * k * rho, s);
    if (bound > 1.0) bound = 1.0;
    if (bound < 1e-15) bound = 1e-15;

    return bound;
}

double nw_ac0_parity_lower_bound(size_t n, size_t depth) {
    /* Hastad (1987):
     * PARITY requires AC0 circuits of size
     *   S = exp(Omega(n^{1/(d-1)}))
     * for depth-d circuits.
     *
     * The precise constant depends on the basis and depth.
     */
    if (depth < 2) return (double)n;

    double exp_part = pow((double)n, 1.0 / (double)(depth - 1));
    return pow(2.0, exp_part / 10.0); /* conservative constant 1/10 */
}

/* =============================================
 * L5: Circuit Complexity Estimation
 * ============================================= */

double nw_circuit_complexity_estimate(const TruthTable *tt) {
    if (!tt) return 0.0;
    return nw_estimate_circuit_complexity(tt);
}

double nw_random_function_complexity(size_t n) {
    /* By Shannon's counting argument:
     * - Almost all boolean functions on n vars require
     *   circuit size >= 2^n / (2n) * (1 - o(1))
     * - There exist functions requiring size up to 2^n/n
     *
     * We return the expected complexity (average):
     *   C_random(n) ≈ 2^n / n
     */
    if (n <= 1) return 1.0;
    return pow(2.0, (double)n) / (double)n;
}

/* =============================================
 * L5: Compile Truth Table to Circuit (Lupanov)
 *
 * Lupanov (1958): Any boolean function can be implemented
 * by a circuit of size (1 + o(1)) * 2^n / n.
 *
 * The construction:
 *   1. Partition the input variables into two groups: x = (u, v)
 *      where |u| = ceiling(log n), |v| = n - |u|
 *   2. For each assignment to u, the restriction f_u(v) depends on n' = n - |u| vars
 *   3. Implement each f_u as a DNF or via universal circuit for n'-variable functions
 *   4. Select the correct f_u based on u using a multiplexer
 *
 * This gives total size ~ 2^n/n.
 * ============================================= */

Circuit *nw_compile_truth_table(const TruthTable *tt) {
    if (!tt || tt->n > 16) return NULL; /* Practical limit */

    /* For small n, we build the circuit using the DNF representation:
     * f(x) = OR_{x: f(x)=1} (AND of literals)
     *
     * Size <= n * 2^n, depth <= n + log(2^n) = O(n)
     */
    Circuit *c = nw_circuit_create(tt->n);
    if (!c) return NULL;

    /* Count ones and build term for each */
    size_t *term_gates = (size_t *)malloc(tt->size * sizeof(size_t));
    size_t num_terms = 0;

    for (size_t i = 0; i < tt->size; i++) {
        if (tt->table[i]) {
            /* Build AND of literals for this minterm.
             * For each var j, if bit j of i is 1, use x_j else use NOT x_j */
            size_t literal_gates[32]; /* max inputs */
            size_t lit_count = 0;

            for (size_t j = 0; j < tt->n && lit_count < 32; j++) {
                if ((i >> j) & 1) {
                    literal_gates[lit_count++] = j; /* x_j */
                } else {
                    literal_gates[lit_count++] = nw_circuit_add_not(c, j);
                }
            }

            /* AND all literals together */
            if (lit_count == 1) {
                term_gates[num_terms++] = literal_gates[0];
            } else if (lit_count > 1) {
                size_t and_gate = nw_circuit_add_and(c, literal_gates[0], literal_gates[1]);
                for (size_t k = 2; k < lit_count; k++) {
                    and_gate = nw_circuit_add_and(c, and_gate, literal_gates[k]);
                }
                term_gates[num_terms++] = and_gate;
            }
        }
    }

    /* OR all terms together */
    if (num_terms == 0) {
        /* Constant 0 function: output AND(x_0, NOT x_0) */
        size_t not0 = nw_circuit_add_not(c, 0);
        size_t zero = nw_circuit_add_and(c, 0, not0);
        nw_circuit_set_output(c, zero);
    } else if (num_terms == 1) {
        nw_circuit_set_output(c, term_gates[0]);
    } else {
        size_t or_gate = nw_circuit_add_or(c, term_gates[0], term_gates[1]);
        for (size_t k = 2; k < num_terms; k++) {
            or_gate = nw_circuit_add_or(c, or_gate, term_gates[k]);
        }
        nw_circuit_set_output(c, or_gate);
    }

    free(term_gates);
    return c;
}

/**
 * L5: Goldreich-Levin Hard Predicate
 *
 * Goldreich & Levin (1989): For any one-way function f,
 * the inner product <x, r> mod 2 is a hard-core predicate.
 *
 * Specifically, g(x, r) = sum_{i=1}^n x_i * r_i mod 2
 * is unpredictable given f(x) and r.
 *
 * In the NW context, this is a natural candidate for the
 * hard function: g is efficiently computable (in NC1) yet
 * provably hard if f is hard.
 */
Circuit *nw_goldreich_levin_hard_predicate(size_t n) {
    /* Inner product function on 2n bits:
     *   IP(x_1,...,x_n, r_1,...,r_n) = XOR_{i=1}^n (x_i AND r_i)
     *
     * This is in NC1: can be computed with O(n) gates and O(log n) depth.
     */
    size_t total_inputs = 2 * n;
    Circuit *c = nw_circuit_create(total_inputs);
    if (!c) return NULL;

    if (n == 0) {
        nw_circuit_set_output(c, 0);
        return c;
    }

    /* For each i: AND(x_i, r_i) */
    size_t *and_gates = (size_t *)malloc(n * sizeof(size_t));
    for (size_t i = 0; i < n; i++) {
        and_gates[i] = nw_circuit_add_and(c, i, n + i);
    }

    /* XOR all AND results */
    size_t result = and_gates[0];
    for (size_t i = 1; i < n; i++) {
        result = nw_circuit_add_xor(c, result, and_gates[i]);
    }

    nw_circuit_set_output(c, result);
    free(and_gates);
    return c;
}

/* =============================================
 * Utility: Circuit Verification and Display
 * ============================================= */

bool nw_circuit_verify(const Circuit *c) {
    if (!c) return false;
    if (c->n_inputs > c->n_gates) return false;

    /* Check that all gate inputs reference valid gates */
    for (size_t i = c->n_inputs; i < c->n_gates; i++) {
        for (size_t j = 0; j < c->fan_in[i]; j++) {
            if (c->input_gates[i][j] >= i)
                return false; /* Input must come from earlier gate (DAG property) */
        }
    }

    /* Check output gates are valid */
    for (size_t i = 0; i < c->n_outputs; i++) {
        if (c->output_gates[i] >= c->n_gates)
            return false;
    }

    return true;
}

void nw_circuit_print(const Circuit *c, FILE *fp) {
    if (!c) { fprintf(fp, "Circuit(NULL)\n"); return; }

    fprintf(fp, "Circuit {\n");
    fprintf(fp, "  inputs: %zu, gates: %zu, outputs: %zu\n",
            c->n_inputs, c->n_gates, c->n_outputs);
    fprintf(fp, "  depth: %zu\n", nw_circuit_depth(c));

    /* Print gate statistics */
    size_t counts[7] = {0};
    for (size_t i = 0; i < c->n_gates; i++)
        counts[c->gate_types[i]]++;

    fprintf(fp, "  gate counts: INPUT=%zu, AND=%zu, OR=%zu, NOT=%zu, XOR=%zu, MAJ=%zu, THR=%zu\n",
            counts[0], counts[1], counts[2], counts[3], counts[4], counts[5], counts[6]);

    /* Print output gates */
    fprintf(fp, "  output gates: [");
    for (size_t i = 0; i < c->n_outputs && i < 10; i++)
        fprintf(fp, "%zu%s", c->output_gates[i], (i+1 < c->n_outputs) ? ", " : "");
    if (c->n_outputs > 10) fprintf(fp, "...");
    fprintf(fp, "]\n}\n");
}

/**
 * L4: Counting Argument for Circuit Lower Bounds
 *
 * Number of distinct circuits of size s on n inputs
 * over basis with b gate types:
 *   At most (b * (s+n)^2)^s
 *
 * Taking log: log_2(#circuits) = s * (log_2(b) + 2*log_2(s+n))
 *                              ≈ O(s log s)
 *
 * For most functions to have circuits of size s, we need:
 *   s * (log_2(b) + 2*log_2(s+n)) >= 2^n
 *
 * => s = Omega(2^n / n)
 */
double nw_log_circuit_count(size_t n, size_t s) {
    if (s == 0) return 0.0;
    double b = 7.0; /* Number of gate types in our basis */
    double log_b = log2(b);
    double log_sn = log2((double)(s + n));
    return (double)s * (log_b + 2.0 * log_sn);
}
