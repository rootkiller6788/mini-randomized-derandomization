/******************************************************************************
 * probabilistic_circuits.h — Probabilistic Circuits & BPP ⊆ P/poly
 *
 * Models randomized computation as Boolean circuits with random inputs,
 * covering:
 *   - Probabilistic circuit families
 *   - Adleman's Theorem: BPP ⊆ P/poly
 *   - Circuit lower bounds for derandomization
 *   - Randomized vs deterministic circuit complexity
 *
 * L1: Definitions — circuit family with random inputs
 * L3: Math Structures — Boolean circuits as DAGs
 * L4: Fundamental Laws — BPP ⊆ P/poly proof
 * L6: Canonical Problems — Circuit-SAT for probabilistic circuits
 *
 * References:
 *   Arora-Barak Ch.6, §7.4
 *   Vollmer, Introduction to Circuit Complexity (1999)
 *   Adleman (1978)
 ******************************************************************************/

#ifndef PROBABILISTIC_CIRCUITS_H
#define PROBABILISTIC_CIRCUITS_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include "randomized_classes.h"

/* ================================================================
 * L1: Core Data Structures
 * ================================================================ */

/** @brief Gate types in a Boolean circuit. */
typedef enum {
    GATE_INPUT   = 0,  /**< Input variable x_i or random bit r_j    */
    GATE_AND     = 1,  /**< AND gate (fan-in 2)                     */
    GATE_OR      = 2,  /**< OR  gate (fan-in 2)                     */
    GATE_NOT     = 3,  /**< NOT gate (fan-in 1)                     */
    GATE_NAND    = 4,  /**< NAND gate (universal)                   */
    GATE_XOR     = 5,  /**< XOR gate                                */
    GATE_MAJ     = 6,  /**< MAJORITY gate (unbounded fan-in)        */
    GATE_THRESHOLD = 7,/**< THRESHOLD_k gate (≥k inputs true)       */
    GATE_OUTPUT  = 8,  /**< Output gate                             */
    GATE_CONST_0 = 9,  /**< Constant 0                              */
    GATE_CONST_1 = 10, /**< Constant 1                              */
} GateType;

/** @brief A single gate in a Boolean circuit. */
typedef struct {
    GateType    type;
    size_t      id;             /**< Unique gate identifier          */
    size_t      fanin;          /**< Number of inputs                */
    size_t     *inputs;         /**< Array of input gate IDs         */
    bool        is_random;      /**< Is this gate a random bit?      */
    size_t      random_index;   /**< If random, which bit index?     */
    bool        is_input_var;   /**< Is this an input variable?      */
    size_t      input_var_idx;  /**< If input var, which x_i?        */
    bool        value;          /**< Computed value (after eval)     */
    bool        evaluated;      /**< Has this gate been evaluated?   */
} CircuitGate;

/** @brief A Boolean circuit with explicit structure. */
typedef struct {
    CircuitGate *gates;         /**< Array of all gates              */
    size_t       num_gates;     /**< Total number of gates           */
    size_t       num_inputs;    /**< Number of input variables x_i   */
    size_t       num_random;    /**< Number of random bit inputs     */
    size_t       num_outputs;   /**< Number of output gates          */
    size_t      *output_gates;  /**< Indices of output gates         */
    size_t       size;          /**< Circuit size = number of gates  */
    size_t       depth;         /**< Longest path from input to output */
} BooleanCircuit;

/** @brief A probabilistic circuit family: {C_n} with random inputs. */
typedef struct {
    BooleanCircuit **circuits;  /**< C_n for each input length        */
    size_t           max_n;     /**< Maximum input length supported   */
    size_t          *sizes;     /**< size(n) for each n              */
    size_t          *random_bits;/**< m(n) random bits per input length */
    size_t          *depths;    /**< depth(n) for each n              */
} ProbabilisticCircuitFamily;

/* ================================================================
 * L2: Circuit Construction API
 * ================================================================ */

/**
 * @brief Create an empty Boolean circuit.
 *
 * @param num_inputs  Number of input variables x_0,...,x_{n-1}
 * @param num_random  Number of random bits r_0,...,r_{m-1}
 * @param num_outputs Number of output bits
 * @param max_gates   Estimated maximum number of gates (for allocation)
 * @return            New circuit (caller must free with circuit_destroy)
 */
BooleanCircuit *circuit_create(
    size_t num_inputs, size_t num_random,
    size_t num_outputs, size_t max_gates);

/**
 * @brief Add a gate to the circuit.
 *
 * @return Gate ID (index into gates array)
 */
size_t circuit_add_gate(BooleanCircuit *c, GateType type,
                        const size_t *inputs, size_t fanin);

/**
 * @brief Add input variable gate.
 */
size_t circuit_add_input(BooleanCircuit *c, size_t var_idx);

/**
 * @brief Add random bit gate.
 */
size_t circuit_add_random_bit(BooleanCircuit *c, size_t bit_idx);

/**
 * @brief Add constant gate.
 */
size_t circuit_add_constant(BooleanCircuit *c, bool value);

/**
 * @brief Set a gate as output.
 */
void circuit_set_output(BooleanCircuit *c, size_t gate_id, size_t output_idx);

/**
 * @brief Evaluate the circuit on given inputs and random bits.
 *
 * @param c          The circuit
 * @param x          Input values (length = num_inputs)
 * @param r          Random bits (length = num_random)
 * @param outputs    Output array (length = num_outputs)
 */
void circuit_evaluate(
    const BooleanCircuit *c,
    const bool *x, const bool *r,
    bool *outputs);

/**
 * @brief Compute circuit depth via topological longest path.
 */
size_t circuit_compute_depth(BooleanCircuit *c);

/**
 * @brief Free circuit resources.
 */
void circuit_destroy(BooleanCircuit *c);

/* ================================================================
 * L3: Circuit Transformations
 * ================================================================ */

/**
 * @brief Convert circuit to use only NAND gates.
 *
 * Any Boolean function can be computed by a circuit of NAND gates
 * with at most polynomial blow-up in size.
 *
 * @return New circuit (caller must free)
 */
BooleanCircuit *circuit_to_nand_basis(const BooleanCircuit *c);

/**
 * @brief Compute circuit size (total number of non-input gates).
 */
size_t circuit_count_gates(const BooleanCircuit *c);

/**
 * @brief Verify circuit is a DAG (no cycles).
 *
 * @return true if acyclic
 */
bool circuit_is_dag(const BooleanCircuit *c);

/**
 * @brief Topological sort of gates for evaluation order.
 *
 * @param order  Output: gate IDs in topological order
 * @return       true if successful (acyclic)
 */
bool circuit_topological_sort(const BooleanCircuit *c, size_t *order);

/* ================================================================
 * L4: Adleman's Theorem — BPP ⊆ P/poly
 * ================================================================ */

/**
 * @brief Demonstrates Adleman's Theorem: given a BPP machine with error
 * ≤ 2^{-2n}, compute the required advice size and verify the construction.
 *
 * Construction: There are at most 2^n inputs of length n. If error per
 * input is ≤ 2^{-2n}, then by union bound, Pr[∃x: M errs on x] ≤ 2^n·2^{-2n}
 * = 2^{-n} < 1. So there exists r that works for ALL inputs.
 *
 * This r is the polynomial advice (size = m(n) = poly(n) random bits).
 *
 * Reference: Arora-Barak Theorem 7.15
 *
 * @param n                Input length
 * @param error_per_input  Error probability per input
 * @param advice_size      Output: m(n) — number of random bits
 * @return                 true if probability of existence > 0
 */
bool adleman_construction(
    size_t n, double error_per_input,
    size_t *advice_size);

/**
 * @brief Simulate Adleman: given a BPP decider and an error bound,
 * find an advice string r (by exhaustive search over random seeds)
 * that works for all inputs in a given test set.
 *
 * @param test_inputs   Array of test input strings
 * @param num_tests     Number of test inputs
 * @param decider       The BPP decider function
 * @param m             Number of random bits
 * @param oracle        True membership oracle
 * @param advice_out    Output: good advice string r
 * @return              true if advice found
 */
bool adleman_find_advice(
    const char **test_inputs, size_t num_tests,
    bool (*decider)(const char*, size_t, RandomSource*),
    size_t m,
    bool (*oracle)(const char*, size_t),
    bool *advice_out);

/* ================================================================
 * L6: Canonical Problems for Probabilistic Circuits
 * ================================================================ */

/**
 * @brief Probabilistic Circuit-SAT:
 * Given a circuit C with n inputs and m random bits, decide whether
 * there exists an input x such that Pr_r[C(x,r)=1] > 1/2.
 *
 * This is a canonical problem for MA (Merlin-Arthur), the
 * probabilistic analog of NP.
 *
 * @param c      The circuit
 * @param trials Number of random trials for probability estimation
 * @return       true if ∃x: acceptance probability > 1/2
 */
bool prob_circuit_sat(const BooleanCircuit *c, size_t trials);

/**
 * @brief Approximate counting version:
 * Estimate the number of x for which Pr[C(x,r)=1] > 1/2.
 *
 * @param c      The circuit
 * @param trials Trials per estimate
 * @param count  Output: estimated count
 * @return       Confidence radius
 */
size_t prob_circuit_count(
    const BooleanCircuit *c, size_t trials,
    double *count);

/* ================================================================
 * L5: Circuit Lower Bound functions
 * ================================================================ */

/**
 * @brief Compute the PARITY function via a circuit.
 *
 * PARITY(x_1,...,x_n) = x_1 ⊕ x_2 ⊕ ... ⊕ x_n
 *
 * PARITY requires exponential size for constant-depth circuits
 * (Håstad's Switching Lemma), but is easy (size O(n)) with depth O(log n).
 *
 * @param n  Number of inputs
 * @return   Circuit computing PARITY
 */
BooleanCircuit *circuit_parity(size_t n);

/**
 * @brief Compute MAJORITY via a circuit (Valiant's construction).
 *
 * MAJORITY(x_1,...,x_n) = 1 iff Σ x_i > n/2
 *
 * MAJORITY has polynomial-size, O(log n) depth circuits.
 */
BooleanCircuit *circuit_majority(size_t n);

/**
 * @brief Probabilistic circuit for Approximate Majority.
 *
 * Given n bits with ≥ 2/3 ones, amplify to ≥ 1 - 2^{-k} agreement.
 * Uses O(n log n) gates.
 */
BooleanCircuit *circuit_approximate_majority(size_t n, size_t k);

/**
 * @brief Compute circuit complexity lower bound via gate-counting argument.
 *
 * (Shannon's counting argument: most Boolean functions require circuits
 * of size Ω(2^n/n).)
 *
 * @param n             Number of inputs
 * @param gate_types    Number of allowed gate types
 * @param max_fanin     Maximum fan-in
 * @return              Lower bound on circuit size for almost all functions
 */
size_t shannon_lower_bound(size_t n, size_t gate_types, size_t max_fanin);

#endif /* PROBABILISTIC_CIRCUITS_H */
