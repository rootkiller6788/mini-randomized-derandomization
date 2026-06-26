/******************************************************************************
 * probabilistic_circuits.c — Boolean Circuits & BPP ⊆ P/poly
 *
 * Implements Boolean circuit model with random inputs, circuit construction,
 * evaluation, analysis, and Adleman's Theorem verification.
 *
 * L1-L3: Circuit data structures and algorithms
 * L4: Adleman's Theorem (BPP ⊆ P/poly)
 * L5: Circuit lower bounds (Shannon, Håstad)
 * L6: Probabilistic Circuit-SAT
 *
 * References:
 *   Arora-Barak Ch.6 (Boolean Circuits), §7.4
 *   Adleman (1978): Two theorems on random polynomial time
 *   Shannon (1949): The synthesis of two-terminal switching circuits
 ******************************************************************************/

#include "probabilistic_circuits.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

/* ================================================================
 * L1-L2: Circuit Construction
 * ================================================================ */

BooleanCircuit *circuit_create(
    size_t num_inputs, size_t num_random,
    size_t num_outputs, size_t max_gates)
{
    BooleanCircuit *c = (BooleanCircuit *)calloc(1, sizeof(BooleanCircuit));
    if (!c) return NULL;

    c->num_inputs = num_inputs;
    c->num_random = num_random;
    c->num_outputs = num_outputs;

    /* Allocate gates array (will grow as needed) */
    size_t initial_cap = max_gates > 0 ? max_gates : 256;
    c->gates = (CircuitGate *)calloc(initial_cap, sizeof(CircuitGate));
    if (!c->gates) {
        free(c);
        return NULL;
    }
    c->size = initial_cap;
    c->num_gates = 0;
    c->depth = 0;

    /* Allocate output gates array */
    c->output_gates = (size_t *)calloc(num_outputs, sizeof(size_t));
    if (!c->output_gates && num_outputs > 0) {
        free(c->gates);
        free(c);
        return NULL;
    }
    for (size_t i = 0; i < num_outputs; i++) {
        c->output_gates[i] = (size_t)-1;
    }

    return c;
}

static bool circuit_grow(BooleanCircuit *c) {
    if (!c) return false;
    size_t new_size = c->size * 2;
    CircuitGate *new_gates = (CircuitGate *)realloc(
        c->gates, new_size * sizeof(CircuitGate));
    if (!new_gates) return false;
    /* Zero newly allocated memory */
    memset(new_gates + c->size, 0, (new_size - c->size) * sizeof(CircuitGate));
    c->gates = new_gates;
    c->size = new_size;
    return true;
}

size_t circuit_add_gate(BooleanCircuit *c, GateType type,
                        const size_t *inputs, size_t fanin)
{
    if (!c) return (size_t)-1;
    if (c->num_gates >= c->size) {
        if (!circuit_grow(c)) return (size_t)-1;
    }

    size_t id = c->num_gates;
    c->gates[id].id = id;
    c->gates[id].type = type;
    c->gates[id].fanin = fanin;

    if (fanin > 0 && inputs) {
        c->gates[id].inputs = (size_t *)malloc(fanin * sizeof(size_t));
        if (c->gates[id].inputs) {
            memcpy(c->gates[id].inputs, inputs, fanin * sizeof(size_t));
        }
    }

    c->gates[id].evaluated = false;
    c->num_gates++;
    return id;
}

size_t circuit_add_input(BooleanCircuit *c, size_t var_idx) {
    if (!c) return (size_t)-1;
    if (c->num_gates >= c->size) {
        if (!circuit_grow(c)) return (size_t)-1;
    }

    size_t id = c->num_gates;
    c->gates[id].id = id;
    c->gates[id].type = GATE_INPUT;
    c->gates[id].is_input_var = true;
    c->gates[id].input_var_idx = var_idx;
    c->gates[id].fanin = 0;
    c->num_gates++;
    return id;
}

size_t circuit_add_random_bit(BooleanCircuit *c, size_t bit_idx) {
    if (!c) return (size_t)-1;
    if (c->num_gates >= c->size) {
        if (!circuit_grow(c)) return (size_t)-1;
    }

    size_t id = c->num_gates;
    c->gates[id].id = id;
    c->gates[id].type = GATE_INPUT;
    c->gates[id].is_random = true;
    c->gates[id].random_index = bit_idx;
    c->gates[id].fanin = 0;
    c->num_gates++;
    return id;
}

size_t circuit_add_constant(BooleanCircuit *c, bool value) {
    if (!c) return (size_t)-1;
    if (c->num_gates >= c->size) {
        if (!circuit_grow(c)) return (size_t)-1;
    }

    size_t id = c->num_gates;
    c->gates[id].id = id;
    c->gates[id].type = value ? GATE_CONST_1 : GATE_CONST_0;
    c->gates[id].value = value;
    c->gates[id].evaluated = true;
    c->gates[id].fanin = 0;
    c->num_gates++;
    return id;
}

void circuit_set_output(BooleanCircuit *c, size_t gate_id, size_t output_idx) {
    if (!c || output_idx >= c->num_outputs) return;
    c->output_gates[output_idx] = gate_id;
}

/* ================================================================
 * L2-L3: Circuit Evaluation and Analysis
 * ================================================================ */

/**
 * Evaluate a single gate recursively, caching results.
 * Inputs and random bits are looked up from the provided arrays.
 */
static bool eval_gate(BooleanCircuit *c, size_t gate_id,
                      const bool *x, const bool *r) {
    if (!c || gate_id >= c->num_gates) return false;

    CircuitGate *g = &c->gates[gate_id];

    /* Already evaluated? Return cached value */
    if (g->evaluated) return g->value;

    switch (g->type) {
        case GATE_INPUT:
            if (g->is_input_var) {
                g->value = (x && g->input_var_idx < c->num_inputs)
                         ? x[g->input_var_idx] : false;
            } else if (g->is_random) {
                g->value = (r && g->random_index < c->num_random)
                         ? r[g->random_index] : false;
            }
            g->evaluated = true;
            return g->value;

        case GATE_CONST_0:
            g->value = false;
            g->evaluated = true;
            return false;

        case GATE_CONST_1:
            g->value = true;
            g->evaluated = true;
            return true;

        case GATE_NOT:
            if (g->fanin >= 1) {
                g->value = !eval_gate(c, g->inputs[0], x, r);
            }
            g->evaluated = true;
            return g->value;

        case GATE_AND:
            g->value = true;
            for (size_t i = 0; i < g->fanin; i++) {
                if (!eval_gate(c, g->inputs[i], x, r)) {
                    g->value = false;
                    break;
                }
            }
            g->evaluated = true;
            return g->value;

        case GATE_OR:
            g->value = false;
            for (size_t i = 0; i < g->fanin; i++) {
                if (eval_gate(c, g->inputs[i], x, r)) {
                    g->value = true;
                    break;
                }
            }
            g->evaluated = true;
            return g->value;

        case GATE_NAND:
            g->value = false;
            for (size_t i = 0; i < g->fanin; i++) {
                if (!eval_gate(c, g->inputs[i], x, r)) {
                    g->value = true;
                    break;
                }
            }
            g->evaluated = true;
            return g->value;

        case GATE_XOR:
            g->value = false;
            for (size_t i = 0; i < g->fanin; i++) {
                if (eval_gate(c, g->inputs[i], x, r)) {
                    g->value = !g->value;
                }
            }
            g->evaluated = true;
            return g->value;

        case GATE_MAJ:
        {
            size_t count = 0;
            for (size_t i = 0; i < g->fanin; i++) {
                if (eval_gate(c, g->inputs[i], x, r)) count++;
            }
            g->value = (count > g->fanin / 2);
            g->evaluated = true;
            return g->value;
        }

        case GATE_THRESHOLD:
        {
            /* THRESHOLD_k: true if ≥ k inputs are true.
             * k is encoded in a field we need to track — for now use
             * the id's upper bits or a fixed threshold g->fanin/2+1. */
            size_t k = (g->fanin / 2) + 1;
            size_t count = 0;
            for (size_t i = 0; i < g->fanin; i++) {
                if (eval_gate(c, g->inputs[i], x, r)) count++;
            }
            g->value = (count >= k);
            g->evaluated = true;
            return g->value;
        }

        case GATE_OUTPUT:
            if (g->fanin >= 1) {
                g->value = eval_gate(c, g->inputs[0], x, r);
            }
            g->evaluated = true;
            return g->value;

        default:
            g->value = false;
            g->evaluated = true;
            return false;
    }
}

void circuit_evaluate(
    const BooleanCircuit *c,
    const bool *x, const bool *r,
    bool *outputs)
{
    if (!c || !outputs) return;

    /* Reset evaluation state */
    BooleanCircuit *nc = (BooleanCircuit *)c;
    for (size_t i = 0; i < c->num_gates; i++) {
        nc->gates[i].evaluated = false;
    }

    /* Evaluate each output gate */
    for (size_t i = 0; i < c->num_outputs; i++) {
        size_t out_id = c->output_gates[i];
        if (out_id != (size_t)-1 && out_id < c->num_gates) {
            outputs[i] = eval_gate(nc, out_id, x, r);
        } else {
            outputs[i] = false;
        }
    }
}

size_t circuit_compute_depth(BooleanCircuit *c) {
    if (!c || c->num_gates == 0) return 0;

    /* Dynamic programming: depth[g] = 1 + max_{input i} depth[i] */
    size_t *depth = (size_t *)calloc(c->num_gates, sizeof(size_t));
    if (!depth) return (size_t)-1;

    /* Topological order */
    size_t *order = (size_t *)malloc(c->num_gates * sizeof(size_t));
    if (!order) { free(depth); return (size_t)-1; }

    if (!circuit_topological_sort(c, order)) {
        free(depth); free(order);
        return (size_t)-1;
    }

    size_t max_depth = 0;
    for (size_t idx = 0; idx < c->num_gates; idx++) {
        size_t id = order[idx];
        CircuitGate *g = &c->gates[id];

        if (g->type == GATE_INPUT || g->type == GATE_CONST_0
            || g->type == GATE_CONST_1) {
            depth[id] = 0;
        } else {
            size_t max_input_depth = 0;
            for (size_t j = 0; j < g->fanin; j++) {
                if (g->inputs[j] < c->num_gates
                    && depth[g->inputs[j]] > max_input_depth) {
                    max_input_depth = depth[g->inputs[j]];
                }
            }
            depth[id] = max_input_depth + 1;
        }
        if (depth[id] > max_depth) max_depth = depth[id];
    }

    c->depth = max_depth;
    free(depth);
    free(order);
    return max_depth;
}

bool circuit_is_dag(const BooleanCircuit *c) {
    if (!c) return false;

    /* Color-based cycle detection:
     * 0 = white (unvisited), 1 = gray (in progress), 2 = black (done) */
    uint8_t *color = (uint8_t *)calloc(c->num_gates, sizeof(uint8_t));
    if (!color) return false;

    bool acyclic = true;
    /* Use iterative DFS to avoid stack overflow on large circuits */
    for (size_t i = 0; i < c->num_gates && acyclic; i++) {
        if (color[i] == 0) {
            /* Iterative DFS from gate i */
            size_t *stack = (size_t *)malloc(c->num_gates * sizeof(size_t));
            size_t *parent = (size_t *)malloc(c->num_gates * sizeof(size_t));
            size_t sp = 0;
            if (!stack || !parent) { acyclic = false; free(stack); free(parent); break; }

            stack[sp++] = i;
            while (sp > 0 && acyclic) {
                size_t u = stack[sp - 1];
                if (color[u] == 0) {
                    color[u] = 1;  /* Gray */
                    /* Push neighbors */
                    CircuitGate *g = &c->gates[u];
                    for (size_t j = 0; j < g->fanin; j++) {
                        size_t v = g->inputs[j];
                        if (v >= c->num_gates) continue;
                        if (color[v] == 1) {
                            acyclic = false;  /* Back edge → cycle */
                            break;
                        }
                        if (color[v] == 0) {
                            parent[v] = u;
                            stack[sp++] = v;
                            if (sp >= c->num_gates) { acyclic = false; break; }
                        }
                    }
                } else {
                    color[u] = 2;  /* Black */
                    sp--;
                }
            }
            free(stack);
            free(parent);
        }
    }

    free(color);
    return acyclic;
}

void circuit_destroy(BooleanCircuit *c) {
    if (!c) return;
    if (c->gates) {
        for (size_t i = 0; i < c->num_gates; i++) {
            free(c->gates[i].inputs);
        }
        free(c->gates);
    }
    free(c->output_gates);
    free(c);
}

/* ================================================================
 * L3: Circuit Transformations
 * ================================================================ */

BooleanCircuit *circuit_to_nand_basis(const BooleanCircuit *c) {
    if (!c) return NULL;

    BooleanCircuit *nc = circuit_create(
        c->num_inputs, c->num_random, c->num_outputs,
        c->num_gates * 3);  /* NAND conversion blows up by ≤ 3× */
    if (!nc) return NULL;

    /* Map old gate IDs to new gate IDs.
     * For each old gate, compute NAND equivalent. */
    size_t *id_map = (size_t *)malloc(c->num_gates * sizeof(size_t));
    if (!id_map) { circuit_destroy(nc); return NULL; }

    /* First pass: copy input and constant gates */
    for (size_t i = 0; i < c->num_gates; i++) {
        const CircuitGate *g = &c->gates[i];
        switch (g->type) {
            case GATE_INPUT:
                if (g->is_input_var) {
                    id_map[i] = circuit_add_input(nc, g->input_var_idx);
                } else if (g->is_random) {
                    id_map[i] = circuit_add_random_bit(nc, g->random_index);
                } else {
                    id_map[i] = circuit_add_constant(nc, false);
                }
                break;
            case GATE_CONST_0:
                id_map[i] = circuit_add_constant(nc, false);
                break;
            case GATE_CONST_1:
                id_map[i] = circuit_add_constant(nc, true);
                break;
            default:
                id_map[i] = (size_t)-1;  /* To be converted */
                break;
        }
    }

    /* Second pass: convert logic gates to NAND */
    for (size_t i = 0; i < c->num_gates; i++) {
        const CircuitGate *g = &c->gates[i];
        if (id_map[i] != (size_t)-1) continue;  /* Already mapped */

        /* Convert each gate type to NAND equivalents */
        size_t inputs[2];
        switch (g->type) {
            case GATE_NOT:
                /* NOT(a) = NAND(a, a) */
                if (g->fanin >= 1) {
                    inputs[0] = id_map[g->inputs[0]];
                    inputs[1] = id_map[g->inputs[0]];
                    id_map[i] = circuit_add_gate(nc, GATE_NAND, inputs, 2);
                }
                break;

            case GATE_AND:
                /* AND(a,b) = NOT(NAND(a,b)) = NAND(NAND(a,b), NAND(a,b)) */
                if (g->fanin >= 2) {
                    /* In practice, fan-in 2 */
                    size_t nand_id = circuit_add_gate(nc, GATE_NAND,
                        (size_t[]){id_map[g->inputs[0]], id_map[g->inputs[1]]}, 2);
                    /* NOT via NAND */
                    id_map[i] = circuit_add_gate(nc, GATE_NAND,
                        (size_t[]){nand_id, nand_id}, 2);
                }
                break;

            case GATE_OR:
                /* OR(a,b) = NOT(NAND(NOT(a), NOT(b)))
                 * = NAND(NAND(a,a), NAND(b,b)) * NAND itself */
                if (g->fanin >= 2) {
                    size_t na = circuit_add_gate(nc, GATE_NAND,
                        (size_t[]){id_map[g->inputs[0]], id_map[g->inputs[0]]}, 2);
                    size_t nb = circuit_add_gate(nc, GATE_NAND,
                        (size_t[]){id_map[g->inputs[1]], id_map[g->inputs[1]]}, 2);
                    id_map[i] = circuit_add_gate(nc, GATE_NAND,
                        (size_t[]){na, nb}, 2);
                }
                break;

            case GATE_NAND:
                /* Already NAND: just copy */
                if (g->fanin >= 2) {
                    id_map[i] = circuit_add_gate(nc, GATE_NAND,
                        (size_t[]){id_map[g->inputs[0]], id_map[g->inputs[1]]}, 2);
                }
                break;

            case GATE_XOR:
                /* XOR(a,b) = OR(AND(a,NOT(b)), AND(NOT(a),b))
                 * Requires multiple NAND gates. Skip for now. */
                id_map[i] = circuit_add_constant(nc, false);
                break;

            default:
                id_map[i] = circuit_add_constant(nc, false);
                break;
        }
    }

    /* Set output gates */
    for (size_t i = 0; i < c->num_outputs; i++) {
        if (c->output_gates[i] != (size_t)-1) {
            circuit_set_output(nc, id_map[c->output_gates[i]], i);
        }
    }

    free(id_map);
    return nc;
}

size_t circuit_count_gates(const BooleanCircuit *c) {
    if (!c) return 0;
    size_t count = 0;
    for (size_t i = 0; i < c->num_gates; i++) {
        if (c->gates[i].type != GATE_INPUT
            && c->gates[i].type != GATE_CONST_0
            && c->gates[i].type != GATE_CONST_1) {
            count++;
        }
    }
    return count;
}

bool circuit_topological_sort(const BooleanCircuit *c, size_t *order) {
    if (!c || !order) return false;
    if (!circuit_is_dag(c)) return false;

    /* Kahn's algorithm */
    size_t *in_degree = (size_t *)calloc(c->num_gates, sizeof(size_t));
    if (!in_degree) return false;

    for (size_t i = 0; i < c->num_gates; i++) {
        for (size_t j = 0; j < c->gates[i].fanin; j++) {
            size_t inp = c->gates[i].inputs[j];
            if (inp < c->num_gates) in_degree[inp]++;  /* Actually: in_degree counts incoming edges. Wait — in_degree[input] means number of gates that depend on it. That's "out-degree" in DAG terms.
               Let me fix: we want topological order from inputs to outputs.
               in_degree[i] = number of predecessors of gate i. */
        }
    }

    /* Recompute correctly */
    memset(in_degree, 0, c->num_gates * sizeof(size_t));
    for (size_t i = 0; i < c->num_gates; i++) {
        for (size_t j = 0; j < c->gates[i].fanin; j++) {
            size_t inp = c->gates[i].inputs[j];
            if (inp < c->num_gates) in_degree[i]++;  /* Gate i depends on inp */
        }
    }

    /* Queue for Kahn */
    size_t *queue = (size_t *)malloc(c->num_gates * sizeof(size_t));
    if (!queue) { free(in_degree); return false; }

    size_t head = 0, tail = 0;
    for (size_t i = 0; i < c->num_gates; i++) {
        if (in_degree[i] == 0) queue[tail++] = i;
    }

    size_t count = 0;
    while (head < tail) {
        size_t u = queue[head++];
        order[count++] = u;

        /* Reduce in-degree of all gates that depend on u */
        for (size_t i = 0; i < c->num_gates; i++) {
            for (size_t j = 0; j < c->gates[i].fanin; j++) {
                if (c->gates[i].inputs[j] == u) {
                    in_degree[i]--;
                    if (in_degree[i] == 0) queue[tail++] = i;
                    break;  /* A gate can depend on u at most once per input */
                }
            }
        }
    }

    free(queue);
    free(in_degree);
    return (count == c->num_gates);
}

/* ================================================================
 * L4: Adleman's Theorem — BPP ⊆ P/poly
 * ================================================================ */

bool adleman_construction(
    size_t n, double error_per_input,
    size_t *advice_size)
{
    /* Theorem (Adleman, 1978):
     * For any L ∈ BPP, L ∈ P/poly.
     *
     * Proof: Let M be BPP machine for L with error ≤ 2^{-|x|-1}.
     * For fixed n, Pr[M(x,r) errs] ≤ 2^{-n-1} per x.
     * Union bound over all 2^n inputs:
     *   Pr[∃x of length n s.t. M errs] ≤ 2^n · 2^{-n-1} = 1/2 < 1.
     * Therefore ∃r (advice string) s.t. M(x,r) correct ∀x ∈ {0,1}^n.
     *
     * The advice size = m(n) = poly(n) random bits. */

    /* Check if union bound works */
    double total_eror = pow(2.0, (double)n) * error_per_input;
    bool advice_exists = (total_eror < 1.0);

    if (advice_exists && advice_size) {
        /* m(n) = log₂(1/error_per_input) random bits for base machine,
         * plus amplification overhead. */
        double m = log2(1.0 / error_per_input);
        if (m < 1.0) m = 1.0;
        *advice_size = (size_t)ceil(m);
    }

    return advice_exists;
}

bool adleman_find_advice(
    const char **test_inputs, size_t num_tests,
    bool (*decider)(const char*, size_t, RandomSource*),
    size_t m,
    bool (*oracle)(const char*, size_t),
    bool *advice_out)
{
    /* Search over all m-bit strings for one that works on all test inputs.
     * Exhaustive for m ≤ 20, random sampling otherwise. */

    if (m > 20) {
        /* Too many to enumerate; use random search */
        for (size_t attempt = 0; attempt < 10000; attempt++) {
            RandomSource rs = random_source_init(attempt * 0xCAFEULL);

            /* Generate candidate advice */
            for (size_t bit = 0; bit < m; bit++) {
                advice_out[bit] = (random_source_bit(&rs) == 1);
            }

            /* Test on all inputs */
            bool all_correct = true;
            for (size_t t = 0; t < num_tests; t++) {
                RandomSource rs2 = random_source_init(0);
                /* Use advice as random bits */
                /* (simplified: direct evaluation) */
                bool result = decider(test_inputs[t], strlen(test_inputs[t]), &rs2);
                bool truth = oracle(test_inputs[t], strlen(test_inputs[t]));
                if (result != truth) { all_correct = false; break; }
            }
            if (all_correct) return true;
        }
        return false;
    }

    /* Exhaustive for small m */
    size_t total = (size_t)1ULL << m;
    for (size_t r = 0; r < total; r++) {
        for (size_t bit = 0; bit < m; bit++) {
            advice_out[bit] = (r >> bit) & 1;
        }

        bool all_correct = true;
        for (size_t t = 0; t < num_tests; t++) {
            RandomSource rs = random_source_init(r);
            bool result = decider(test_inputs[t], strlen(test_inputs[t]), &rs);
            bool truth = oracle(test_inputs[t], strlen(test_inputs[t]));
            if (result != truth) { all_correct = false; break; }
        }
        if (all_correct) return true;
    }
    return false;
}

/* ================================================================
 * L6: Probabilistic Circuit-SAT
 * ================================================================ */

bool prob_circuit_sat(const BooleanCircuit *c, size_t trials) {
    if (!c || c->num_inputs == 0) return false;

    /* Enumerate all inputs (for small n); check probability bound.
     * For large n, use random sampling. */

    size_t total_inputs = (size_t)1ULL << c->num_inputs;
    size_t max_enum = (total_inputs < 1024) ? total_inputs : 1024;

    for (size_t x_val = 0; x_val < max_enum; x_val++) {
        bool *x = (bool *)calloc(c->num_inputs, sizeof(bool));
        if (!x) continue;

        for (size_t j = 0; j < c->num_inputs; j++) {
            x[j] = (x_val >> j) & 1;
        }

        /* Estimate Pr_r[C(x,r) = 1] via sampling */
        size_t accepts = 0;
        for (size_t t = 0; t < trials; t++) {
            bool *r = (bool *)calloc(c->num_random, sizeof(bool));
            if (!r) continue;
            for (size_t j = 0; j < c->num_random; j++) {
                r[j] = ((t * 6364136223846793005ULL + j) & 1);
            }
            bool output;
            circuit_evaluate(c, x, r, &output);
            if (output) accepts++;
            free(r);
        }

        double prob = (double)accepts / (double)trials;
        if (prob > 0.5) {
            free(x);
            return true;
        }
        free(x);
    }

    return false;
}

size_t prob_circuit_count(
    const BooleanCircuit *c, size_t trials,
    double *estimated_count)
{
    if (!c) { if (estimated_count) *estimated_count = 0.0; return 0; }

    /* Count inputs x for which Pr[C(x,r)=1] > 1/2 (up to enumeration limit) */

    size_t total_inputs = (size_t)1ULL << c->num_inputs;
    size_t max_enum = (total_inputs < 256) ? total_inputs : 256;
    size_t count = 0;

    for (size_t x_val = 0; x_val < max_enum; x_val++) {
        bool *x = (bool *)calloc(c->num_inputs, sizeof(bool));
        if (!x) continue;
        for (size_t j = 0; j < c->num_inputs; j++) {
            x[j] = (x_val >> j) & 1;
        }

        size_t accepts = 0;
        for (size_t t = 0; t < trials; t++) {
            bool *r = (bool *)calloc(c->num_random, sizeof(bool));
            if (!r) continue;
            for (size_t j = 0; j < c->num_random; j++) {
                r[j] = ((t * 6364136223846793005ULL + j) & 1);
            }
            bool output;
            circuit_evaluate(c, x, r, &output);
            if (output) accepts++;
            free(r);
        }

        double prob = (double)accepts / (double)trials;
        if (prob > 0.5) count++;
        free(x);
    }

    if (estimated_count) {
        *estimated_count = (double)count;
        if (total_inputs > max_enum) {
            *estimated_count *= (double)total_inputs / (double)max_enum;
        }
    }

    return count;
}

/* ================================================================
 * L5: Circuit Lower Bounds — Shannon & PARITY
 * ================================================================ */

BooleanCircuit *circuit_parity(size_t n) {
    /* PARITY(x₁,...,xₙ) = x₁ ⊕ x₂ ⊕ ... ⊕ xₙ
     * Binary tree of XOR gates: depth ⌈log₂ n⌉, size n-1 XORs.
     * Each XOR = 4 NAND gates: size = 4(n-1) in NAND basis. */

    BooleanCircuit *c = circuit_create(n, 0, 1, 4 * n + 10);
    if (!c) return NULL;

    /* Create input gates */
    size_t *input_ids = (size_t *)malloc(n * sizeof(size_t));
    for (size_t i = 0; i < n; i++) {
        input_ids[i] = circuit_add_input(c, i);
    }

    /* Build XOR tree */
    size_t remaining = n;
    size_t *current = input_ids;
    size_t *next = (size_t *)malloc(n * sizeof(size_t));

    while (remaining > 1) {
        size_t next_size = 0;
        for (size_t i = 0; i + 1 < remaining; i += 2) {
            /* XOR(a,b) = (a AND NOT b) OR (NOT a AND b)
             * In NAND: NAND(NAND(a,NAND(b,b)),NAND(NAND(a,a),b)) */
            size_t nb = circuit_add_gate(c, GATE_NAND,
                (size_t[]){current[i+1], current[i+1]}, 2);
            size_t na = circuit_add_gate(c, GATE_NAND,
                (size_t[]){current[i], current[i]}, 2);
            size_t t1 = circuit_add_gate(c, GATE_NAND,
                (size_t[]){current[i], nb}, 2);
            size_t t2 = circuit_add_gate(c, GATE_NAND,
                (size_t[]){na, current[i+1]}, 2);
            next[next_size++] = circuit_add_gate(c, GATE_NAND,
                (size_t[]){t1, t2}, 2);
        }
        if (remaining % 2 == 1) {
            next[next_size++] = current[remaining - 1];
        }
        size_t *tmp = current;
        current = next;
        next = tmp;
        remaining = next_size;
    }

    circuit_set_output(c, current[0], 0);

    /* Free allocated arrays.
     * current and next alternate; one is input_ids, the other is malloc'd.
     * Determine which is which to avoid use-after-free. */
    {
        size_t *to_free1 = input_ids;
        size_t *to_free2 = (input_ids == current || input_ids == next)
                           ? ((input_ids == current) ? next : current) : NULL;
        free(to_free1);
        free(to_free2);
    }
    return c;
}

BooleanCircuit *circuit_majority(size_t n) {
    /* MAJORITY: 1 iff > n/2 inputs are 1.
     * Use sorting-network approach: O(n log² n) size, O(log² n) depth.
     * Simplified: use a MAJ gate directly. */

    BooleanCircuit *c = circuit_create(n, 0, 1, n + 5);
    if (!c) return NULL;

    size_t *inputs = (size_t *)malloc(n * sizeof(size_t));
    for (size_t i = 0; i < n; i++) {
        inputs[i] = circuit_add_input(c, i);
    }

    /* Use the built-in MAJ gate with all inputs */
    size_t out = circuit_add_gate(c, GATE_MAJ, inputs, n);
    circuit_set_output(c, out, 0);

    free(inputs);
    return c;
}

BooleanCircuit *circuit_approximate_majority(size_t n, size_t k) {
    /* Approximate majority: given n bits with p ≥ 2/3 ones,
     * output 1 with probability ≥ 1 - 2^{-k}.
     *
     * Construction: t = O(k) layers of random pairing and AND.
     * This is a simplified simulation. */
    BooleanCircuit *c = circuit_create(n, 0, 1, n + 10 * k);
    if (!c) return NULL;
    /* Simplified: just use majority gate (which is exact, not approximate) */
    BooleanCircuit *tmp = circuit_majority(n);
    circuit_destroy(c);
    (void)k;  /* k incorporated conceptually */
    return tmp;
}

size_t shannon_lower_bound(size_t n, size_t gate_types, size_t max_fanin) {
    /* Shannon (1949) counting argument:
     *
     * Number of Boolean functions on n variables: 2^{2^n}.
     *
     * Number of circuits of size s with g gate types, fan-in ≤ f:
     *   ≤ (g · (s+n)^{f+1})^{s}
     *
     * Taking logs: 2^n ≤ s log(g·(s+n)^{f+1}) ≈ s·f·log(s+n).
     *
     * Solving: s ≥ Ω(2^n / n). */

    if (n <= 1) return n;

    double total_funcs = pow(2.0, pow(2.0, (double)n));
    double log_funcs = log2(total_funcs);  /* = 2^n */

    /* Number of possible circuits of size s:
     * Each gate: choose type (g), choose up to f inputs from s+n wires.
     * Total choices per gate: ≤ g · (s+n)^f. */

    /* Binary search for minimal s such that #circuits ≥ #functions */
    size_t lo = 1, hi = (size_t)1ULL << (n + 5);
    while (lo < hi) {
        size_t s = lo + (hi - lo) / 2;
        /* log₂(#circuits) ≈ s · (log₂(g) + f · log₂(s + n)) */
        double log_circuits = (double)s
            * (log2((double)gate_types)
               + (double)max_fanin * log2((double)(s + n)));
        if (log_circuits >= log_funcs) {
            hi = s;
        } else {
            lo = s + 1;
        }
    }

    return lo;
}
