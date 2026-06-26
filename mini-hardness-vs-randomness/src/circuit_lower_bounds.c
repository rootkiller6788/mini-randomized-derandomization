/******************************************************************************
 * circuit_lower_bounds.c �� Circuit Lower Bounds Implementation
 *
 * Implements Boolean circuit construction, evaluation, and lower bound
 * computation for the hardness-vs-randomness paradigm.
 *
 * Knowledge: L1 (circuit, gate definitions), L3 (DAG structure, evaluation),
 * L4 (Shannon, H?stad, Razborov, Smolensky theorems)
 *
 * References:
 *   Shannon (1949), H?stad (1986), Razborov (1985), Smolensky (1987),
 *   Furst-Saxe-Sipser (1984), Yao (1985), Alon-Boppana (1987)
 *   Arora & Barak (2009), Chapters 6 and 14
 ******************************************************************************/

#include "circuit_lower_bounds.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

/* ================================================================
 * L3: Boolean Circuit Construction and Evaluation
 * ================================================================ */

BooleanCircuit *clb_circuit_create(size_t n) {
    BooleanCircuit *c = malloc(sizeof(BooleanCircuit));
    if (!c) return NULL;
    c->n = n;
    c->num_gates = 0;
    c->num_outputs = 0;
    /* Start with space for input variables + modest working room */
    c->gates = malloc((n + 64) * sizeof(Gate));
    if (!c->gates) { free(c); return NULL; }
    c->output_ids = malloc(16 * sizeof(size_t));
    if (!c->output_ids) { free(c->gates); free(c); return NULL; }
    c->depth = 0;
    c->size = 0;
    /* Create input variable gates */
    for (size_t i = 0; i < n; i++) {
        c->gates[i].type = GATE_INPUT;
        c->gates[i].input1 = i;
        c->gates[i].input2 = i;
        c->gates[i].negated = false;
    }
    c->num_gates = n;
    return c;
}

void clb_circuit_free(BooleanCircuit *circuit) {
    if (!circuit) return;
    free(circuit->gates);
    free(circuit->output_ids);
    free(circuit);
}

size_t clb_circuit_add_gate(BooleanCircuit *circuit, GateType type,
                             size_t in1, size_t in2, bool negated) {
    if (!circuit) return SIZE_MAX;
    /* For unary gates, only in1 matters */
    if (type == GATE_NOT || type == GATE_INPUT) {
        in2 = in1;
    }
    size_t id = circuit->num_gates;
    /* Reallocate if needed */
    size_t cap = (id < 64) ? circuit->n + 64 : id * 2;
    if (id >= cap) {
        size_t new_cap = id * 2;
        Gate *new_gates = realloc(circuit->gates, new_cap * sizeof(Gate));
        if (!new_gates) return SIZE_MAX;
        circuit->gates = new_gates;
    }
    circuit->gates[id].type = type;
    circuit->gates[id].input1 = in1;
    circuit->gates[id].input2 = in2;
    circuit->gates[id].negated = negated;
    circuit->num_gates = id + 1;
    circuit->size = id + 1 - circuit->n; /* size counts non-input gates */
    return id;
}

size_t clb_circuit_set_output(BooleanCircuit *circuit, size_t gate_id) {
    if (!circuit || gate_id >= circuit->num_gates) return circuit->num_outputs;
    size_t idx = circuit->num_outputs;
    /* Reallocate output_ids if needed */
    if (idx % 16 == 0) {
        size_t *new_ids = realloc(circuit->output_ids, (idx + 16) * sizeof(size_t));
        if (!new_ids) return idx;
        circuit->output_ids = new_ids;
    }
    circuit->output_ids[idx] = gate_id;
    circuit->num_outputs = idx + 1;
    return circuit->num_outputs;
}

bool clb_circuit_evaluate(const BooleanCircuit *circuit,
                           const bool *x, bool *result) {
    if (!circuit || !x || !result) return false;
    size_t n = circuit->n;
    size_t N = circuit->num_gates;
    /* Allocate value array on stack for small circuits, heap for large */
    bool *values;
    bool stack_alloc = (N <= 4096);
    if (stack_alloc) {
        values = alloca(N * sizeof(bool));
    } else {
        values = malloc(N * sizeof(bool));
        if (!values) return false;
    }
    /* Topological evaluation: gates are stored in sorted order */
    for (size_t i = 0; i < N; i++) {
        bool v;
        switch (circuit->gates[i].type) {
            case GATE_INPUT: {
                size_t var = circuit->gates[i].input1;
                v = (var < n) ? x[var] : false;
                break;
            }
            case GATE_NOT: {
                size_t in = circuit->gates[i].input1;
                v = (in < i) ? !values[in] : false;
                break;
            }
            case GATE_AND: {
                size_t a = circuit->gates[i].input1;
                size_t b = circuit->gates[i].input2;
                v = (a < i && b < i) ? (values[a] && values[b]) : false;
                break;
            }
            case GATE_OR: {
                size_t a = circuit->gates[i].input1;
                size_t b = circuit->gates[i].input2;
                v = (a < i && b < i) ? (values[a] || values[b]) : false;
                break;
            }
            case GATE_XOR: {
                size_t a = circuit->gates[i].input1;
                size_t b = circuit->gates[i].input2;
                v = (a < i && b < i) ? (values[a] != values[b]) : false;
                break;
            }
            case GATE_MAJORITY:
            case GATE_MODp:
            default:
                v = false;
                break;
        }
        if (circuit->gates[i].negated) v = !v;
        values[i] = v;
    }
    /* Write outputs */
    for (size_t o = 0; o < circuit->num_outputs; o++) {
        size_t gate_id = circuit->output_ids[o];
        result[o] = (gate_id < N) ? values[gate_id] : false;
    }
    if (!stack_alloc) free(values);
    return true;
}

size_t clb_circuit_compute_depth(BooleanCircuit *circuit) {
    if (!circuit) return 0;
    size_t N = circuit->num_gates;
    size_t *depths = calloc(N, sizeof(size_t));
    if (!depths) return 0;
    /* Inputs have depth 0 */
    for (size_t i = 0; i < circuit->n; i++) {
        depths[i] = 0;
    }
    size_t max_depth = 0;
    for (size_t i = circuit->n; i < N; i++) {
        size_t d1 = (circuit->gates[i].input1 < i) ? depths[circuit->gates[i].input1] : 0;
        size_t d2 = (circuit->gates[i].input2 < i) ? depths[circuit->gates[i].input2] : 0;
        size_t d = (d1 > d2 ? d1 : d2) + 1;
        depths[i] = d;
        if (d > max_depth) max_depth = d;
    }
    free(depths);
    circuit->depth = max_depth;
    return max_depth;
}

bool clb_circuit_verify_structure(const BooleanCircuit *circuit) {
    if (!circuit) return false;
    if (circuit->n == 0) return false;
    if (circuit->num_gates < circuit->n) return false;
    /* Every non-input gate must reference earlier gates */
    for (size_t i = circuit->n; i < circuit->num_gates; i++) {
        if (circuit->gates[i].input1 >= i) return false;
        if (circuit->gates[i].type != GATE_NOT && circuit->gates[i].type != GATE_INPUT) {
            if (circuit->gates[i].input2 >= i) return false;
        }
    }
    for (size_t o = 0; o < circuit->num_outputs; o++) {
        if (circuit->output_ids[o] >= circuit->num_gates) return false;
    }
    return true;
}

/* ================================================================
 * L4: Shannon's Counting Lower Bound
 *
 * Theorem (Shannon, 1949): For n sufficiently large, there exists
 * a Boolean function f: {0,1}^n -> {0,1} requiring circuits of
 * size at least 2^n / (3n).
 *
 * Proof via counting: Number of circuits of size �� S is �� ((n+S)^2)^S.
 * Number of Boolean functions is 2^{2^n}.
 * If ((n+S)^2)^S < 2^{2^n}, some function requires size > S.
 * Taking log: S * 2log(n+S) < 2^n => S < 2^n/(2log(n+S)).
 * Asymptotically: S ~ 2^n/n.
 * ================================================================ */

CircuitLB clb_shannon_counting(size_t n) {
    CircuitLB lb;
    lb.n = n;
    lb.theorem = "Shannon (1949) Counting Argument";
    lb.problem = "Most n-variable Boolean functions";
    if (n < 2) {
        lb.gates = 1;
        lb.lower_bound = 1.0;
        return lb;
    }
    /* Bound: S �� 2^n / (3n) for most functions */
    double two_n = pow(2.0, (double)n);
    double bound = two_n / (3.0 * (double)n);
    lb.gates = (size_t)bound;
    lb.lower_bound = bound;
    return lb;
}

/*
 * Lupanov (1958) upper bound:
 * Every Boolean function on n variables can be computed by a circuit
 * of size (1 + o(1)) * 2^n / n.
 * This shows Shannon's lower bound is asymptotically tight.
 */
CircuitLB clb_lupanov_upper(size_t n) {
    CircuitLB lb;
    lb.n = n;
    lb.theorem = "Lupanov (1958) Upper Bound";
    lb.problem = "All n-variable Boolean functions";
    if (n < 2) {
        lb.gates = 1;
        lb.lower_bound = 1.0;
        return lb;
    }
    double two_n = pow(2.0, (double)n);
    double bound = two_n / (double)n;
    /* Lupanov: (1 + o(1)) * 2^n/n �� 1.01 * 2^n/n for large n */
    if (n < 10) {
        bound *= 1.5; /* small n correction */
    }
    lb.gates = (size_t)bound;
    lb.lower_bound = bound;
    return lb;
}

/* ================================================================
 * H?stad's Switching Lemma (1986):
 *
 * PARITY requires AC0 circuits of size �� 2^{c * n^{1/(d-1)}} for
 * depth-d circuits over AND/OR/NOT.
 *
 * Key technique: Random restriction + switching lemma.
 * A depth-d circuit becomes depth-(d-1) with high probability
 * under random restriction, while PARITY remains essentially PARITY.
 * ================================================================ */

CircuitLB clb_hastad_switching(size_t n, size_t depth) {
    CircuitLB lb;
    lb.n = n;
    lb.theorem = "Hastad (1986) Switching Lemma";
    lb.problem = "PARITY not in AC0";
    if (depth < 2) depth = 2;
    if (n < 2) {
        lb.gates = 1;
        lb.lower_bound = 1.0;
        return lb;
    }
    /* Bound: 2^{c * n^{1/(d-1)}} with c �� 0.05 */
    double exponent = pow((double)n, 1.0 / (double)(depth - 1));
    double bound = pow(2.0, 0.05 * exponent) / 10.0;
    lb.gates = (size_t)bound;
    lb.lower_bound = bound;
    return lb;
}

/* ================================================================
 * Smolensky's Lower Bound (1987):
 *
 * MODp (parity mod p) requires exponential-size AC0[q] circuits
 * when p and q are distinct primes.
 *
 * Uses polynomial approximation over GF(q):
 * AC0[q] circuits compute low-degree polynomials over GF(q),
 * but MODp requires degree ��(��n) over GF(q) when p �� q.
 * ================================================================ */

CircuitLB clb_smolensky_modp(size_t n, size_t p) {
    (void)p; /* p parameter reserved for distinguishing MODp from MODq */
    CircuitLB lb;
    lb.n = n;
    lb.theorem = "Smolensky (1987) Algebraic Lower Bound";
    lb.problem = "MODp not in AC0[q]";
    if (n < 2) {
        lb.gates = 1;
        lb.lower_bound = 1.0;
        return lb;
    }
    /* Bound: 2^{c * sqrt{n}} for some constant c �� 0.1 */
    double bound = pow(2.0, 0.1 * sqrt((double)n));
    lb.gates = (size_t)bound;
    lb.lower_bound = bound;
    return lb;
}

/* ================================================================
 * Razborov's Monotone Circuit Lower Bound (1985):
 *
 * CLIQUE_k_n requires monotone circuits of size 2^{��(k^{1/3})}.
 *
 * Method: Method of approximations (sunflower lemma).
 * Construct an "approximator" for CLIQUE that captures the
 * structure while being much smaller than any monotone circuit.
 * ================================================================ */

CircuitLB clb_razborov_clique(size_t n) {
    CircuitLB lb;
    lb.n = n;
    lb.theorem = "Razborov (1985) Method of Approximations";
    lb.problem = "CLIQUE requires large monotone circuits";
    if (n < 3) {
        lb.gates = 1;
        lb.lower_bound = 1.0;
        return lb;
    }
    /* Take k = n^{1/4} for simplicity */
    double k = pow((double)n, 0.25);
    /* Bound: 2^{c * k^{1/3}} with c �� 0.03 */
    double bound = pow(2.0, 0.03 * pow(k, 1.0/3.0));
    lb.gates = (size_t)bound;
    lb.lower_bound = bound;
    return lb;
}

/*
 * Alon-Boppana (1987) refinement of Razborov's bound.
 * Tighter sunflower constructions give better constants.
 */
CircuitLB clb_alon_boppana(size_t n, size_t clique_size) {
    CircuitLB lb;
    lb.n = n;
    lb.theorem = "Alon-Boppana (1987) Refined Monotone Bound";
    lb.problem = "CLIQUE monotone circuit lower bound";
    if (n < 3 || clique_size < 2 || clique_size > n) {
        lb.gates = 1;
        lb.lower_bound = 1.0;
        return lb;
    }
    double k = (double)clique_size;
    double bound = pow(2.0, 0.04 * pow(k, 1.0/3.0));
    lb.gates = (size_t)bound;
    lb.lower_bound = bound;
    return lb;
}

/* ================================================================
 * Majority Function Lower Bound (Furst-Saxe-Sipser 1984, Yao 1985):
 *
 * MAJORITY: {0,1}^n -> {0,1}, MAJ(x) = 1 iff sum(x_i) > n/2.
 * MAJORITY requires AC0 circuits of size exp(n^{1/(d-1)}).
 * ================================================================ */

CircuitLB clb_majority(size_t n) {
    CircuitLB lb;
    lb.n = n;
    lb.theorem = "Furst-Saxe-Sipser (1984) / Yao (1985)";
    lb.problem = "MAJORITY not in AC0";
    if (n < 2) {
        lb.gates = 1;
        lb.lower_bound = 1.0;
        return lb;
    }
    /* For any constant depth d, size �� 2^{n^{1/(d-1)}/poly(n)} */
    double bound = pow(2.0, 0.1 * pow((double)n, 0.33) * (double)n / 10.0);
    if (bound < 2.0) bound = 2.0;
    lb.gates = (size_t)bound;
    lb.lower_bound = bound;
    return lb;
}

/*
 * PARITY circuit lower bound (AC0 version).
 * PARITY(x) = XOR of all bits.
 * PARITY ? AC0: requires super-polynomial size for constant depth.
 */
CircuitLB clb_parity(size_t n) {
    CircuitLB lb;
    lb.n = n;
    lb.theorem = "Hastad (1986) Parity Lower Bound";
    lb.problem = "PARITY not in AC0";
    if (n < 2) {
        lb.gates = 1;
        lb.lower_bound = 1.0;
        return lb;
    }
    /* For constant depth d, size �� 2^{n^{1/(d-1)}/20} */
    double bound = pow(2.0, 0.1 * pow((double)n, 0.5));
    if (bound < 2.0) bound = 2.0;
    lb.gates = (size_t)bound;
    lb.lower_bound = bound;
    return lb;
}

/*
 * Formula lower bound (Andreev/Hastad).
 * An explicit function requiring formulas of size ��(n^{3}).
 * The function: f(g, x) = PARITY(g(x)) where g encodes a selector.
 */
CircuitLB clb_formula_lower_bound(size_t n) {
    CircuitLB lb;
    lb.n = n;
    lb.theorem = "Andreev (1987) / Hastad (1998)";
    lb.problem = "Explicit function with cubic formula size";
    if (n < 4) {
        lb.gates = 1;
        lb.lower_bound = 1.0;
        return lb;
    }
    double bound = pow((double)n, 3.0) / 1000.0;
    lb.gates = (size_t)bound;
    lb.lower_bound = bound;
    return lb;
}

/*
 * Diagonalization lower bound:
 * There exists a function requiring circuits of size > S
 * (non-constructive, via counting/diagonalization).
 */
CircuitLB clb_diagonalization_bound(size_t n, size_t S) {
    CircuitLB lb;
    lb.n = n;
    lb.theorem = "Diagonalization (Complexity Theory)";
    lb.problem = "Existence of hard functions";
    /* The bound S itself is the claim: there exists f with
       circuit complexity > S */
    lb.gates = S;
    lb.lower_bound = (double)S;
    return lb;
}

/* ================================================================
 * L2: Lower Bound Utility Functions
 * ================================================================ */

bool clb_is_strong_enough(const CircuitLB *clb, HRLevel level) {
    if (!clb) return false;
    double lb = clb->lower_bound;
    double n = (double)clb->n;
    switch (level) {
        case HR_WEAK:
            /* Sub-exponential: n^{log n} */
            return lb > pow(n, log(n));
        case HR_MODERATE:
            /* 2^{n^��} for �� > 0.1 */
            return lb > pow(2.0, pow(n, 0.1));
        case HR_STRONG:
            /* 2^{c*n} for c > 0.1 */
            return lb > pow(2.0, 0.1 * n);
        case HR_EXPONENTIAL:
            /* 2^{(1-��)n} for �� < 0.5 */
            return lb > pow(2.0, 0.5 * n);
        default:
            return false;
    }
}

double clb_quality_metric(const CircuitLB *clb) {
    if (!clb || clb->n == 0) return 0.0;
    double n = (double)clb->n;
    /* Quality = log2(lower_bound) / n */
    double log_lb = (clb->lower_bound > 1.0) ? log2(clb->lower_bound) : 0.0;
    return log_lb / n;
}

const char *clb_class_name(CircuitClass cls) {
    switch (cls) {
        case CLB_AC0:      return "AC0";
        case CLB_ACC0:     return "ACC0";
        case CLB_TC0:      return "TC0";
        case CLB_NC1:      return "NC1";
        case CLB_AC1:      return "AC1";
        case CLB_NC:       return "NC";
        case CLB_P_POLY:   return "P/poly";
        case CLB_MONOTONE: return "Monotone";
        case CLB_EXP:      return "EXP-size";
        default:           return "Unknown";
    }
}

const char *clb_type_name(const CircuitLB *clb) {
    if (!clb) return "NULL";
    return clb->theorem ? clb->theorem : "Unspecified";
}

const char *clb_level_name(HRLevel level) {
    switch (level) {
        case HR_WEAK:        return "WEAK";
        case HR_MODERATE:    return "MODERATE";
        case HR_STRONG:      return "STRONG";
        case HR_EXPONENTIAL: return "EXPONENTIAL";
        default:             return "UNKNOWN";
    }
}

size_t clb_min_gates_for_hardness(double hardness, size_t n) {
    /* hardness parameter = log(lower_bound) / n, so
       lower_bound = 2^{hardness * n} */
    double bound = pow(2.0, hardness * (double)n);
    return (size_t)bound;
}

double clb_min_lower_bound_for_level(HRLevel level, size_t n) {
    double dn = (double)n;
    switch (level) {
        case HR_WEAK:
            return pow(dn, log(dn));
        case HR_MODERATE:
            return pow(2.0, pow(dn, 0.1));
        case HR_STRONG:
            return pow(2.0, 0.1 * dn);
        case HR_EXPONENTIAL:
            return pow(2.0, 0.5 * dn);
        default:
            return 1.0;
    }
}

/* ================================================================
 * L3: Concrete Hard Function Construction
 * ================================================================ */

bool *clb_build_parity_truth_table(size_t n) {
    if (n > 24) return NULL; /* Too large: 2^24 = 16M entries */
    size_t rows = (size_t)1 << n;
    bool *tt = malloc(rows * sizeof(bool));
    if (!tt) return NULL;
    for (size_t i = 0; i < rows; i++) {
        /* Population count parity */
        size_t x = i;
        x ^= x >> 1;
        x ^= x >> 2;
        x = (x & 0x1111111111111111ULL) * 0x1111111111111111ULL;
        tt[i] = ((x >> 60) & 1) ? true : false;
    }
    return tt;
}

bool *clb_build_majority_truth_table(size_t n) {
    if (n > 24) return NULL;
    size_t rows = (size_t)1 << n;
    bool *tt = malloc(rows * sizeof(bool));
    if (!tt) return NULL;
    for (size_t i = 0; i < rows; i++) {
        /* Count set bits */
        size_t cnt = 0;
        for (size_t b = 0; b < n; b++) {
            if (i & ((size_t)1 << b)) cnt++;
        }
        tt[i] = (cnt > n / 2);
    }
    return tt;
}

bool *clb_build_random_truth_table(size_t n, uint64_t seed) {
    if (n > 24) return NULL;
    size_t rows = (size_t)1 << n;
    bool *tt = malloc(rows * sizeof(bool));
    if (!tt) return NULL;
    /* xorshift64* PRNG */
    uint64_t state = (seed == 0) ? 0xDEADBEEFCAFEBABEULL : seed;
    for (size_t i = 0; i < rows; i++) {
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
        uint64_t r = state * 0x2545F4914F6CDD1DULL;
        tt[i] = (r & 1) ? true : false;
    }
    return tt;
}

bool clb_heuristic_hard_function_test(const bool *tt, size_t n, size_t S) {
    (void)S; /* reserved for circuit size threshold parameter */
    if (!tt || n > 20) return false;
    /* Heuristic: count distinct patterns of size log(S) and
       check if the function appears complex.
       Crude test: measure entropy of k-bit windows. */
    if (n < 4) return false;
    size_t rows = (size_t)1 << n;
    /* Check if the truth table has balanced output */
    size_t ones = 0;
    for (size_t i = 0; i < rows; i++) {
        if (tt[i]) ones++;
    }
    double balance = (double)ones / (double)rows;
    /* Functions extremely close to constant are easy */
    if (balance < 0.01 || balance > 0.99) return false;
    /* Check autocorrelation with shift-1 */
    size_t matches = 0;
    for (size_t i = 0; i < rows - 1; i++) {
        if (tt[i] == tt[i+1]) matches++;
    }
    double corr = (double)matches / (double)(rows - 1);
    /* High correlation suggests simplicity */
    if (corr > 0.9 || corr < 0.1) return false;
    /* Crude estimate: function passes basic hardness heuristic */
    return true;
}

double clb_estimate_circuit_complexity(const bool *tt, size_t n) {
    if (!tt || n > 20 || n == 0) return 0.0;
    size_t rows = (size_t)1 << n;
    /* Count number of truth-table entries that differ from linear
       threshold functions. Higher discrepancy = higher complexity.
       This is a heuristic, not a rigorous lower bound. */
    size_t disc = 0;
    /* Test against parity */
    size_t parity_diff = 0;
    for (size_t i = 0; i < rows; i++) {
        bool parity = false;
        for (size_t b = 0; b < n; b++) {
            if (i & ((size_t)1 << b)) parity = !parity;
        }
        if (tt[i] != parity) parity_diff++;
    }
    /* Test against monotone threshold */
    size_t thresh_diff = 0;
    for (size_t i = 0; i < rows; i++) {
        size_t cnt = 0;
        for (size_t b = 0; b < n; b++) {
            if (i & ((size_t)1 << b)) cnt++;
        }
        if (tt[i] != (cnt > n/2)) thresh_diff++;
    }
    /* Minimum discrepancy against simple function families gives heuristic */
    disc = (parity_diff < thresh_diff) ? parity_diff : thresh_diff;
    double disc_ratio = (double)disc / (double)rows;
    /* Map to approximate circuit size: higher discrepancy = larger circuits */
    double est = disc_ratio * pow(2.0, (double)n) / (double)n;
    return est;
}

/*
 * Diagonal hard function:
 * Define f(x) = 1 if the x-th circuit of size �� S outputs 0 on input x,
 *             = 0 otherwise.
 * This is a variant of the standard diagonalization proof that
 * there exist functions with high circuit complexity.
 *
 * Since we can't actually enumerate all circuits (there are too many),
 * we approximate by a simpler diagonal construction:
 * f(x) = parity of bits of x at positions determined by S.
 */
bool *clb_build_diagonal_hard_function(size_t n, size_t S) {
    if (n > 24) return NULL;
    size_t rows = (size_t)1 << n;
    bool *tt = malloc(rows * sizeof(bool));
    if (!tt) return NULL;
    /* Use a deterministic construction that breaks low-complexity patterns */
    uint64_t state = S;
    for (size_t i = 0; i < rows; i++) {
        /* Combine bits using S as the mixing parameter */
        uint64_t x = i;
        x ^= state;
        x ^= x >> 17;
        x ^= x << 31;
        x ^= x >> 8;
        x *= 0x9E3779B97F4A7C15ULL;
        tt[i] = (x >> 63) & 1;
    }
    return tt;
}

/*
 * Nisan-Wigderson hard function template:
 * Take a hard predicate g: {0,1}^k -> {0,1} with circuit complexity S,
 * embed into n-variable space via a combinatorial design.
 * The embedded function has circuit complexity related to S.
 */
bool *clb_build_nw_hard_function(size_t n, size_t k, const bool *hard_pred_tt) {
    if (!hard_pred_tt || n > 24 || k > n) return NULL;
    size_t rows = (size_t)1 << n;
    size_t pred_size = (size_t)1 << k;
    bool *tt = malloc(rows * sizeof(bool));
    if (!tt) return NULL;
    for (size_t i = 0; i < rows; i++) {
        /* Extract k bits from i to index into hard predicate.
           Use a simple design: take bits [0..k-1]. */
        size_t pred_idx = i & (pred_size - 1);
        tt[i] = hard_pred_tt[pred_idx];
    }
    return tt;
}