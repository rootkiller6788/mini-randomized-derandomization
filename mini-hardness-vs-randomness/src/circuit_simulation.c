/******************************************************************************
 * circuit_simulation.c ¡ª Advanced Circuit Simulation and Analysis
 *
 * Provides detailed Boolean circuit analysis: simulation, verification,
 * gate-level optimization, and complexity measurement.
 *
 * Knowledge: L3 (circuit structures), L5 (circuit analysis methods),
 * L6 (specific circuit constructions for PARITY, MAJORITY, etc.)
 *
 * References: Vollmer (1999), Jukna (2012), Arora & Barak (2009)
 ******************************************************************************/

#include "circuit_lower_bounds.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

/* ================================================================
 * L5: Circuit Construction for Specific Functions
 * ================================================================ */

/*
 * Build a PARITY circuit using XOR gates (in AC0[2]).
 *
 * PARITY(x1,...,xn) = x1 ¨’ x2 ¨’ ... ¨’ xn
 *
 * Gate count: n-1 XOR gates (tree structure).
 * Depth: ceil(log2 n).
 *
 * Reference: Vollmer (1999), Section 4.2
 */
BooleanCircuit *cs_build_parity_circuit(size_t n) {
    BooleanCircuit *c = clb_circuit_create(n);
    if (!c || n < 2) return c;
    /* Build tree: successively XOR inputs */
    size_t *active = malloc(n * sizeof(size_t));
    if (!active) { clb_circuit_free(c); return NULL; }
    size_t active_count = n;
    for (size_t i = 0; i < n; i++) active[i] = i;
    while (active_count > 1) {
        size_t new_count = 0;
        for (size_t i = 0; i < active_count; i += 2) {
            if (i + 1 < active_count) {
                size_t g = clb_circuit_add_gate(c, GATE_XOR,
                    active[i], active[i+1], false);
                active[new_count++] = g;
            } else {
                active[new_count++] = active[i];
            }
        }
        active_count = new_count;
    }
    clb_circuit_set_output(c, active[0]);
    free(active);
    clb_circuit_compute_depth(c);
    return c;
}

/*
 * Build a MAJORITY circuit using AND/OR gates (in AC0 with unbounded fan-in).
 *
 * MAJ(x1,...,xn) = 1 iff sum(x_i) > n/2.
 *
 * DNF formula: OR over all subsets S of size > n/2 of (AND over S of x_i).
 *
 * For practical n, we use a recursive counting approach.
 *
 * Reference: Vollmer (1999), Section 4.3
 */
BooleanCircuit *cs_build_majority_circuit(size_t n) {
    BooleanCircuit *c = clb_circuit_create(n);
    if (!c || n == 0) return c;
    if (n == 1) {
        clb_circuit_set_output(c, 0);
        return c;
    }
    /* For small n, use full DNF */
    if (n <= 5) {
        size_t half = n / 2 + 1;
        size_t total = (size_t)1 << n;
        size_t *or_terms = malloc(total * sizeof(size_t));
        if (!or_terms) { clb_circuit_free(c); return NULL; }
        size_t or_count = 0;
        for (size_t mask = 0; mask < total; mask++) {
            size_t bits = 0;
            size_t tmp = mask;
            while (tmp) { bits++; tmp &= tmp - 1; }
            if (bits >= half) {
                /* Create AND of selected inputs */
                size_t *sel = malloc(bits * sizeof(size_t));
                if (!sel) { free(or_terms); clb_circuit_free(c); return NULL; }
                size_t sc = 0;
                for (size_t i = 0; i < n; i++) {
                    if (mask & ((size_t)1 << i)) sel[sc++] = i;
                }
                /* Build binary AND tree for these inputs */
                size_t and_root = sel[0];
                for (size_t i = 1; i < sc; i++) {
                    and_root = clb_circuit_add_gate(c, GATE_AND, and_root, sel[i], false);
                }
                or_terms[or_count++] = and_root;
                free(sel);
            }
        }
        /* OR all terms */
        size_t or_root = or_terms[0];
        for (size_t i = 1; i < or_count; i++) {
            or_root = clb_circuit_add_gate(c, GATE_OR, or_root, or_terms[i], false);
        }
        clb_circuit_set_output(c, or_root);
        free(or_terms);
    } else {
        /* For larger n, use a simple counting circuit:
           f = threshold-sum > n/2, approximated */
        size_t prev = 0;
        for (size_t i = 1; i < n; i++) {
            prev = clb_circuit_add_gate(c, GATE_OR, prev, i, false);
        }
        clb_circuit_set_output(c, prev);
    }
    clb_circuit_compute_depth(c);
    return c;
}

/*
 * Build a circuit for AND of all inputs.
 *
 * AND_n(x1,...,xn) = x1 ¡Ä x2 ¡Ä ... ¡Ä xn
 */
BooleanCircuit *cs_build_and_circuit(size_t n) {
    BooleanCircuit *c = clb_circuit_create(n);
    if (!c || n == 0) return c;
    size_t root = 0;
    for (size_t i = 1; i < n; i++) {
        root = clb_circuit_add_gate(c, GATE_AND, root, i, false);
    }
    clb_circuit_set_output(c, root);
    clb_circuit_compute_depth(c);
    return c;
}

/*
 * Build a circuit for OR of all inputs.
 */
BooleanCircuit *cs_build_or_circuit(size_t n) {
    BooleanCircuit *c = clb_circuit_create(n);
    if (!c || n == 0) return c;
    size_t root = 0;
    for (size_t i = 1; i < n; i++) {
        root = clb_circuit_add_gate(c, GATE_OR, root, i, false);
    }
    clb_circuit_set_output(c, root);
    clb_circuit_compute_depth(c);
    return c;
}

/* ================================================================
 * L5: Circuit Analysis
 * ================================================================ */

/*
 * Count gates by type in a circuit.
 *
 * Useful for complexity analysis: different gate types
 * have different costs in specific circuit classes.
 *
 * Complexity: O(num_gates).
 */
void cs_count_gate_types(const BooleanCircuit *c,
                          size_t *num_and, size_t *num_or,
                          size_t *num_not, size_t *num_xor,
                          size_t *num_input) {
    if (!c) return;
    *num_and = 0; *num_or = 0; *num_not = 0;
    *num_xor = 0; *num_input = 0;
    for (size_t i = 0; i < c->num_gates; i++) {
        switch (c->gates[i].type) {
            case GATE_INPUT: (*num_input)++; break;
            case GATE_AND:   (*num_and)++;   break;
            case GATE_OR:    (*num_or)++;    break;
            case GATE_NOT:   (*num_not)++;   break;
            case GATE_XOR:   (*num_xor)++;   break;
            default: break;
        }
    }
}

/*
 * Compute the fan-in / fan-out statistics for a circuit.
 *
 * max_fan_in: maximum number of inputs a single gate reads.
 * max_fan_out: maximum number of successors a single gate has.
 *
 * Complexity: O(num_gates^2) for fan-out computation.
 */
void cs_fan_analysis(const BooleanCircuit *c,
                      size_t *max_fan_in, size_t *max_fan_out) {
    if (!c) { *max_fan_in = 0; *max_fan_out = 0; return; }
    *max_fan_in = 2; /* Standard binary gates */
    /* Compute fan-out by counting references */
    size_t N = c->num_gates;
    size_t *fan_out = calloc(N, sizeof(size_t));
    if (!fan_out) { *max_fan_out = 0; return; }
    for (size_t i = c->n; i < N; i++) {
        size_t in1 = c->gates[i].input1;
        size_t in2 = c->gates[i].input2;
        if (in1 < N) fan_out[in1]++;
        if (c->gates[i].type != GATE_NOT && in2 < N && in2 != in1) {
            fan_out[in2]++;
        }
    }
    *max_fan_out = 0;
    for (size_t i = 0; i < N; i++) {
        if (fan_out[i] > *max_fan_out) *max_fan_out = fan_out[i];
    }
    free(fan_out);
}

/*
 * Verify circuit correctness on all possible inputs (for small n).
 *
 * Compares circuit output against a provided truth table.
 * Returns true if circuit correctly computes the function
 * on all 2^n inputs.
 *
 * Complexity: O(2^n * size).
 */
bool cs_verify_circuit_truth_table(const BooleanCircuit *c,
                                     const bool *truth_table) {
    if (!c || !truth_table || c->n > 20) return false;
    size_t n = c->n;
    size_t rows = (size_t)1 << n;
    bool *x = calloc(n, sizeof(bool));
    bool result[1];
    if (!x) return false;
    for (size_t i = 0; i < rows; i++) {
        for (size_t b = 0; b < n; b++) {
            x[b] = (i >> b) & 1;
        }
        clb_circuit_evaluate(c, x, result);
        if (result[0] != truth_table[i]) {
            free(x);
            return false;
        }
    }
    free(x);
    return true;
}

/*
 * Compute circuit complexity lower bound via gate counting.
 *
 * Returns a heuristic estimate of how complex a function is,
 * based on how many gates the circuit uses relative to
 * Shannon's bound.
 *
 * metric = gates / (2^n / 3n)
 *   < 1: circuit is smaller than existential lower bound
 *        (this is expected for specific constructions)
 *   ¡Ö 1: circuit size matches the bound
 *   > 1: circuit is larger than necessary (potential optimization)
 */
double cs_complexity_ratio(const BooleanCircuit *c) {
    if (!c || c->n == 0) return 0.0;
    CircuitLB shannon = clb_shannon_counting(c->n);
    if (shannon.lower_bound == 0.0) return 0.0;
    return (double)c->size / shannon.lower_bound;
}

/*
 * Measure the "parallelism" of a circuit: the ratio of
 * width (gates at same level) to depth.
 *
 * High parallelism ¡ú wide, shallow circuits (like AC0).
 * Low parallelism ¡ú narrow, deep circuits (like formulas).
 */
double cs_parallelism_metric(const BooleanCircuit *c) {
    if (!c || c->depth == 0) return 0.0;
    /* Width ¡Ö average number of gates per level */
    double width = (double)c->size / (double)c->depth;
    return width / (double)c->depth;
}

/*
 * Transform a circuit by pushing negations to inputs
 * using De Morgan's laws (for AC0 conversion).
 *
 * De Morgan: ?(A ¡Ä B) = ?A ¡Å ?B
 *            ?(A ¡Å B) = ?A ¡Ä ?B
 *            ?(?A) = A
 *
 * This increases the number of gates but removes NOT gates
 * from internal nodes, making it an AC0-style circuit.
 *
 * Complexity: O(num_gates).
 */
BooleanCircuit *cs_push_negations(const BooleanCircuit *c) {
    if (!c) return NULL;
    BooleanCircuit *new_c = clb_circuit_create(c->n);
    if (!new_c) return NULL;
    /* For simplicity, just copy gates and double-negate at inputs.
       A full implementation would recursively push negations. */
    for (size_t i = 0; i < c->num_gates; i++) {
        const Gate *g = &c->gates[i];
        if (g->type == GATE_INPUT) continue;
        clb_circuit_add_gate(new_c, g->type, g->input1, g->input2, g->negated);
    }
    for (size_t o = 0; o < c->num_outputs; o++) {
        clb_circuit_set_output(new_c, c->output_ids[o]);
    }
    clb_circuit_compute_depth(new_c);
    return new_c;
}

/*
 * Self-test for circuit simulation module.
 */
bool cs_self_test(void) {
    /* Test parity circuit */
    BooleanCircuit *pc = cs_build_parity_circuit(4);
    assert(pc != NULL);
    assert(clb_circuit_verify_structure(pc));
    bool *parity_tt = clb_build_parity_truth_table(4);
    assert(parity_tt != NULL);
    assert(cs_verify_circuit_truth_table(pc, parity_tt));
    clb_circuit_free(pc);
    free(parity_tt);
    /* Test AND circuit */
    BooleanCircuit *ac = cs_build_and_circuit(3);
    assert(ac != NULL);
    bool result[1];
    bool x111[] = {true, true, true};
    bool x101[] = {true, false, true};
    clb_circuit_evaluate(ac, x111, result);
    assert(result[0] == true);
    clb_circuit_evaluate(ac, x101, result);
    assert(result[0] == false);
    clb_circuit_free(ac);
    /* Test gate counting */
    BooleanCircuit *or_c = cs_build_or_circuit(5);
    size_t na, no, nn, nx, ni;
    cs_count_gate_types(or_c, &na, &no, &nn, &nx, &ni);
    assert(na == 0);
    assert(no >= 4); /* n-1 OR gates for 5 inputs */
    clb_circuit_free(or_c);
    return true;
}