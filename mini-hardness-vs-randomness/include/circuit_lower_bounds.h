#ifndef CLB_H
#define CLB_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>

/* Gate types for Boolean circuits */
typedef enum {
    GATE_INPUT = 0,
    GATE_AND,
    GATE_OR,
    GATE_NOT,
    GATE_XOR,
    GATE_MAJORITY,
    GATE_MODp
} GateType;

/* A single gate in a Boolean circuit DAG */
typedef struct {
    GateType type;
    size_t input1;
    size_t input2;
    bool negated;
} Gate;

/* Boolean circuit: DAG of gates computing f: {0,1}^n -> {0,1}^m */
typedef struct {
    size_t n;
    size_t num_gates;
    size_t num_outputs;
    Gate *gates;
    size_t *output_ids;
    size_t depth;
    size_t size;
} BooleanCircuit;

/* Circuit lower bound certificate */
typedef struct {
    size_t n;
    size_t gates;
    double lower_bound;
    const char *theorem;
    const char *problem;
} CircuitLB;

/* Complexity class for circuit families */
typedef enum {
    CLB_AC0,
    CLB_ACC0,
    CLB_TC0,
    CLB_NC1,
    CLB_AC1,
    CLB_NC,
    CLB_P_POLY,
    CLB_MONOTONE,
    CLB_EXP
} CircuitClass;

/* Hardness level derived from circuit lower bounds */
typedef enum {
    HR_WEAK = 0,
    HR_MODERATE,
    HR_STRONG,
    HR_EXPONENTIAL
} HRLevel;

/* Circuit construction and evaluation */
BooleanCircuit *clb_circuit_create(size_t n);
void clb_circuit_free(BooleanCircuit *circuit);
size_t clb_circuit_add_gate(BooleanCircuit *circuit, GateType type, size_t in1, size_t in2, bool negated);
size_t clb_circuit_set_output(BooleanCircuit *circuit, size_t gate_id);
bool clb_circuit_evaluate(const BooleanCircuit *circuit, const bool *x, bool *result);
size_t clb_circuit_compute_depth(BooleanCircuit *circuit);
bool clb_circuit_verify_structure(const BooleanCircuit *circuit);

/* L4: Circuit Lower Bound Theorems */
CircuitLB clb_shannon_counting(size_t n);
CircuitLB clb_lupanov_upper(size_t n);
CircuitLB clb_hastad_switching(size_t n, size_t depth);
CircuitLB clb_smolensky_modp(size_t n, size_t p);
CircuitLB clb_razborov_clique(size_t n);
CircuitLB clb_alon_boppana(size_t n, size_t clique_size);
CircuitLB clb_majority(size_t n);
CircuitLB clb_parity(size_t n);
CircuitLB clb_formula_lower_bound(size_t n);
CircuitLB clb_diagonalization_bound(size_t n, size_t S);

/* Utility functions */
bool clb_is_strong_enough(const CircuitLB *clb, HRLevel level);
double clb_quality_metric(const CircuitLB *clb);
const char *clb_class_name(CircuitClass cls);
const char *clb_type_name(const CircuitLB *clb);
const char *clb_level_name(HRLevel level);
size_t clb_min_gates_for_hardness(double hardness, size_t n);
double clb_min_lower_bound_for_level(HRLevel level, size_t n);

/* Concrete hard function construction */
bool *clb_build_parity_truth_table(size_t n);
bool *clb_build_majority_truth_table(size_t n);
bool *clb_build_random_truth_table(size_t n, uint64_t seed);
bool clb_heuristic_hard_function_test(const bool *tt, size_t n, size_t S);
double clb_estimate_circuit_complexity(const bool *tt, size_t n);
bool *clb_build_diagonal_hard_function(size_t n, size_t S);
bool *clb_build_nw_hard_function(size_t n, size_t k, const bool *hard_pred_tt);

#endif /* CLB_H */