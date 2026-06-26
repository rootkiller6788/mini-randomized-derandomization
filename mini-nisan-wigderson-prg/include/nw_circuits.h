/**
 * @file nw_circuits.h
 * @brief Boolean circuit model for the NW PRG derandomization framework.
 *
 * This module provides:
 *   1. Boolean circuit representation and evaluation
 *   2. Circuit lower bound constructions (PARITY, MAJORITY)
 *   3. Circuit complexity analysis
 *   4. Håstad's switching lemma for AC0 lower bounds
 *   5. Circuit-based distinguisher construction
 *
 * The circuit model is fundamental to the NW construction because
 * hardness is defined with respect to circuit size. The smaller the
 * circuit that can compute a function, the "easier" it is.
 *
 * Reference:
 *   Arora & Barak, Ch 6 (Boolean circuits)
 *   Vollmer, "Introduction to Circuit Complexity", 1999
 *   Håstad, "Computational Limitations of Small-Depth Circuits", 1987
 *   Razborov, "Lower bounds for monotone circuit complexity", 1985
 */

#ifndef NW_CIRCUITS_H
#define NW_CIRCUITS_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include "nw_core.h"

/* ──────────────────────────────────────────────
 * L1: Circuit Gate Types
 * ────────────────────────────────────────────── */

typedef enum {
    NW_GATE_INPUT = 0,      /**< Input gate (variable or constant) */
    NW_GATE_AND   = 1,      /**< AND gate (fan-in 2) */
    NW_GATE_OR    = 2,      /**< OR gate (fan-in 2) */
    NW_GATE_NOT   = 3,      /**< NOT gate (fan-in 1) */
    NW_GATE_XOR   = 4,      /**< XOR gate (fan-in 2) */
    NW_GATE_MAJ   = 5,      /**< MAJORITY gate (unbounded fan-in) */
    NW_GATE_THR   = 6       /**< Threshold gate (unbounded fan-in) */
} NWGateType;

/* ──────────────────────────────────────────────
 * L1/L3: Circuit Structures
 * ────────────────────────────────────────────── */

/**
 * @brief Full-featured boolean circuit with explicit DAG structure.
 *
 * L1: Definition — Boolean circuit
 * L3: Mathematical Structure — Directed acyclic graph
 */
typedef struct Circuit {
    size_t n_inputs;        /**< Number of primary inputs */
    size_t n_gates;         /**< Total number of gates (including inputs) */
    size_t n_outputs;       /**< Number of outputs */
    size_t n_wires;         /**< Total number of wires */
    NWGateType *gate_types; /**< Type of each gate */
    size_t *fan_in;         /**< Number of inputs to each gate */
    size_t **input_gates;   /**< input_gates[g][i] = index of i-th input */
    size_t *output_gates;   /**< Which gates produce outputs */
    bool *values;           /**< Evaluation workspace */
} Circuit;

/* ──────────────────────────────────────────────
 * L5: Circuit Construction
 * ────────────────────────────────────────────── */

/**
 * @brief Create an empty circuit.
 *
 * @param n_inputs  Number of primary inputs
 * @return Newly allocated circuit
 */
Circuit *nw_circuit_create(size_t n_inputs);

/**
 * @brief Add an AND gate to the circuit.
 * @param c  Circuit
 * @param in1  First input gate index
 * @param in2  Second input gate index
 * @return Index of the new gate
 */
size_t nw_circuit_add_and(Circuit *c, size_t in1, size_t in2);

/**
 * @brief Add an OR gate to the circuit.
 * @param c  Circuit
 * @param in1  First input gate index
 * @param in2  Second input gate index
 * @return Index of the new gate
 */
size_t nw_circuit_add_or(Circuit *c, size_t in1, size_t in2);

/**
 * @brief Add a NOT gate to the circuit.
 * @param c  Circuit
 * @param in1  Input gate index
 * @return Index of the new gate
 */
size_t nw_circuit_add_not(Circuit *c, size_t in1);

/**
 * @brief Add an XOR gate using AND/OR/NOT decomposition.
 *
 * x ⊕ y = (x ∧ ¬y) ∨ (¬x ∧ y)
 *
 * @param c  Circuit
 * @param in1  First input
 * @param in2  Second input
 * @return Index of the new XOR output gate
 */
size_t nw_circuit_add_xor(Circuit *c, size_t in1, size_t in2);

/**
 * @brief Add a MAJORITY gate with unbounded fan-in.
 *
 * MAJ(x_1, ..., x_k) = 1 iff ∑ x_i ≥ k/2.
 *
 * Implemented by AKS sorting-network-based construction:
 * O(k log k) AND/OR/NOT gates.
 *
 * @param c       Circuit
 * @param inputs  Array of input gate indices
 * @param k       Number of inputs
 * @return Index of the MAJ output gate
 */
size_t nw_circuit_add_majority(Circuit *c, const size_t *inputs, size_t k);

/**
 * @brief Set a gate as output.
 * @param c       Circuit
 * @param gate_idx  Gate to designate as output
 */
void nw_circuit_set_output(Circuit *c, size_t gate_idx);

/**
 * @brief Free a circuit.
 * @param c  Circuit to free
 */
void nw_circuit_free(Circuit *c);

/* ──────────────────────────────────────────────
 * L5: Circuit Evaluation
 * ────────────────────────────────────────────── */

/**
 * @brief Evaluate a circuit on a given input.
 *
 * Processes gates in topological order.
 * Complexity: O(n_gates)
 *
 * @param c      Circuit
 * @param input  n_inputs boolean values
 * @param output Output buffer (size n_outputs)
 */
void nw_circuit_evaluate(const Circuit *c, const bool *input, bool *output);

/**
 * @brief Evaluate a single-output circuit.
 * @param c      Circuit (must have n_outputs = 1)
 * @param input  Input vector
 * @return The single output value
 */
bool nw_circuit_eval_single(const Circuit *c, const bool *input);

/**
 * @brief Count circuit size.
 * @param c  Circuit
 * @return Total number of gates (excluding inputs)
 */
size_t nw_circuit_gate_count(const Circuit *c);

/**
 * @brief Compute circuit depth.
 *
 * Depth of input gate = 0.
 * Depth of gate = 1 + max(depth of inputs).
 *
 * @param c  Circuit
 * @return Maximum gate depth
 */
size_t nw_circuit_depth(const Circuit *c);

/* ──────────────────────────────────────────────
 * L6: Canonical Functions with Circuit Constructions
 * ────────────────────────────────────────────── */

/**
 * @brief Build a circuit computing PARITY of n bits.
 *
 * PARITY(x) = x_1 ⊕ x_2 ⊕ ... ⊕ x_n.
 *
 * Naive: O(n) gates using XOR tree.
 * But PARITY ∉ AC0 (Ajtai 1983, Furst-Saxe-Sipser 1984).
 * In the basis {AND, OR, NOT}, PARITY requires size 2^{n^{Ω(1/d)}}
 * for depth-d circuits (Håstad 1987).
 *
 * L6: Canonical Problem — PARITY function
 * L4: Fundamental Law — PARITY ∉ AC0
 *
 * @param n  Number of inputs
 * @return Circuit computing PARITY
 */
Circuit *nw_circuit_parity(size_t n);

/**
 * @brief Build a circuit computing MAJORITY of n bits.
 *
 * MAJORITY(x) = 1 iff ∑ x_i ≥ n/2.
 *
 * MAJORITY ∈ TC0 (can be done with threshold gates).
 * In AC0, MAJORITY requires size 2^{n^{Ω(1/d)}} depth-d circuits.
 *
 * L6: Canonical Problem — MAJORITY function
 *
 * @param n  Number of inputs
 * @return Circuit computing MAJORITY
 */
Circuit *nw_circuit_majority(size_t n);

/**
 * @brief Build a circuit computing the MOD-p function.
 *
 * MOD-p(x) = 1 iff ∑ x_i ≡ 0 (mod p).
 *
 * Razborov-Smolensky: MOD-p ∉ AC0[p'] for distinct primes p, p'.
 *
 * L6: Canonical Problem — MOD-p function
 * L8: Advanced Topic — Razborov-Smolensky lower bounds
 *
 * @param n  Number of inputs
 * @param p  Modulus
 * @return Circuit computing MOD-p
 */
Circuit *nw_circuit_mod_p(size_t n, size_t p);

/* ──────────────────────────────────────────────
 * L4: Circuit Lower Bounds
 * ────────────────────────────────────────────── */

/**
 * @brief Shannon's counting lower bound.
 *
 * At most 2^{2s log s + O(s)} distinct circuits of size s exist.
 * There are 2^{2^n} boolean functions on n variables.
 * Therefore most functions require size ≥ 2^n/(2n).
 *
 * L4: Fundamental Law — Shannon's lower bound (1949)
 *
 * @param n  Number of inputs
 * @return Lower bound on circuit size for almost all functions
 */
size_t nw_shannon_circuit_lower_bound(size_t n);

/**
 * @brief Lupanov's upper bound: every function has circuits
 * of size (1+o(1))·2^n/n.
 *
 * @param n  Number of inputs
 * @return Upper bound
 */
double nw_lupanov_upper_bound_val(size_t n);

/* ──────────────────────────────────────────────
 * L8: Håstad Switching Lemma
 * ────────────────────────────────────────────── */

/**
 * @brief Check if Håstad's switching lemma applies to given parameters.
 *
 * Håstad 1987: A k-DNF (k-CNF) under random restriction with
 * parameter ρ becomes a decision tree of small depth with high
 * probability, which can be expressed as a CNF (DNF).
 *
 * L8: Advanced Topic — Håstad's Switching Lemma
 *
 * @param depth  Circuit depth
 * @param n      Number of inputs
 * @return true if switching lemma provides non-trivial bound
 */
bool nw_hastad_applies(size_t depth, size_t n);

/**
 * @brief Compute the switching lemma probability bound.
 *
 * Pr[a k-DNF cannot be expressed as an (s-1)-CNF after ρ-random restriction]
 * ≤ (α_k)^s where α_k is a constant depending on k and ρ.
 *
 * @param n       Number of variables
 * @param depth   Circuit depth
 * @param rho     1 - ρ = fraction of variables left unfixed
 * @return Switching probability bound
 */
double nw_switching_lemma_bound(size_t n, size_t depth, double rho);

/**
 * @brief Compute AC0 lower bound from switching lemma.
 *
 * PARITY requires AC0 circuits of size exp(n^{Θ(1/(d-1))})
 * for depth-d circuits.
 *
 * @param n      Number of inputs
 * @param depth  Circuit depth
 * @return size lower bound
 */
double nw_ac0_parity_lower_bound(size_t n, size_t depth);

/* ──────────────────────────────────────────────
 * L5: Circuit Complexity Analysis
 * ────────────────────────────────────────────── */

/**
 * @brief Estimate circuit complexity from truth table.
 *
 * Uses the mass production formula: most functions require
 * approximately 2^n/n gates.
 *
 * @param tt  Truth table
 * @return Estimated minimum circuit size
 */
double nw_circuit_complexity_estimate(const TruthTable *tt);

/**
 * @brief Compute the circuit complexity of a random function.
 *
 * By Shannon's counting argument, a random function has
 * circuit complexity Θ(2^n/n).
 *
 * @param n  Number of inputs
 * @return Expected circuit complexity
 */
double nw_random_function_complexity(size_t n);

/* ──────────────────────────────────────────────
 * L5: Compile Hard Function into Circuit
 * ────────────────────────────────────────────── */

/**
 * @brief Compile a hard function (as truth table) into a circuit
 * using the Lupanov decomposition.
 *
 * This builds the asymptotically optimal circuit of size ~ 2^n/n.
 *
 * @param tt  Truth table
 * @return Circuit computing the function
 */
Circuit *nw_compile_truth_table(const TruthTable *tt);

/**
 * @brief Build a circuit that serves as a candidate hard function
 * for the NW construction.
 *
 * Uses the Goldreich-Levin hard predicate construction:
 * f(x, r) = ⟨x, r⟩ mod 2 (inner product) which is known to be
 * a hard-core bit of any one-way function.
 *
 * @param n  Input length (half of total; total = 2n)
 * @return Circuit computing the inner product
 */
Circuit *nw_goldreich_levin_hard_predicate(size_t n);

/* ──────────────────────────────────────────────
 * Utility: Circuit Verification
 * ────────────────────────────────────────────── */

/**
 * @brief Verify that a circuit is syntactically valid.
 *
 * Checks: no cycles, gate indices in range, fan-in correct.
 *
 * @param c  Circuit to verify
 * @return true if valid
 */
bool nw_circuit_verify(const Circuit *c);

/**
 * @brief Print a circuit in a human-readable format.
 * @param c   Circuit
 * @param fp  Output stream
 */
void nw_circuit_print(const Circuit *c, FILE *fp);

/**
 * @brief Count the number of circuits of size s on n inputs.
 *
 * Used for counting argument analysis.
 *
 * @param n  Number of inputs
 * @param s  Circuit size
 * @return log2(number of circuits of size s on n inputs)
 */
double nw_log_circuit_count(size_t n, size_t s);

#endif /* NW_CIRCUITS_H */
