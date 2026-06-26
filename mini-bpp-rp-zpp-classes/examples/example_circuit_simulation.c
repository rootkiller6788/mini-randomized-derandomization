/******************************************************************************
 * example_circuit_simulation.c — Probabilistic Circuit Family Demo
 *
 * Demonstrates Adleman's Theorem (BPP ⊆ P/poly) by constructing
 * a probabilistic circuit family and showing how it can be
 * derandomized with polynomial advice.
 *
 * L4: Fundamental Laws — BPP ⊆ P/poly (Adleman, 1978)
 * L6: Canonical Problems — Probabilistic Circuit-SAT
 *
 * References: Adleman (1978), Arora-Barak Theorem 7.15
 ******************************************************************************/

#include "probabilistic_circuits.h"
#include "randomized_classes.h"
#include "chernoff_bounds.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void) {
    printf("=== Probabilistic Circuits & Adleman's Theorem Demo ===\n\n");

    /* Part 1: Construct and evaluate a simple probabilistic circuit */
    printf("--- Part 1: Boolean Circuit Construction ---\n\n");

    /* C(x₁,x₂; r₁) = (x₁ AND r₁) OR (x₂ AND NOT r₁)
     * This circuit computes a randomized function.
     * With r=0: outputs x₂; with r=1: outputs x₁ */
    BooleanCircuit *c = circuit_create(2, 1, 1, 20);

    size_t x0 = circuit_add_input(c, 0);
    size_t x1 = circuit_add_input(c, 1);
    size_t r0 = circuit_add_random_bit(c, 0);

    /* NOT r0 */
    size_t not_r0 = circuit_add_gate(c, GATE_NOT, (size_t[]){r0}, 1);

    /* t1 = x1 AND r0 */
    size_t t1 = circuit_add_gate(c, GATE_AND, (size_t[]){x0, r0}, 2);

    /* t2 = x2 AND NOT r0 */
    size_t t2 = circuit_add_gate(c, GATE_AND, (size_t[]){x1, not_r0}, 2);

    /* out = t1 OR t2 */
    size_t out = circuit_add_gate(c, GATE_OR, (size_t[]){t1, t2}, 2);

    circuit_set_output(c, out, 0);

    printf("Circuit C(x₀,x₁; r₀) = (x₀ ∧ r₀) ∨ (x₁ ∧ ¬r₀)\n");
    printf("Gates: %zu  Size: %zu\n\n", c->num_gates, circuit_count_gates(c));

    printf("Truth table:\n");
    printf(" x₀ x₁ | r₀=0 | r₀=1\n");
    printf("-------+------+------\n");

    for (int v = 0; v < 4; v++) {
        bool x[2] = {(v>>0)&1, (v>>1)&1};
        bool r[2] = {false};
        bool out0, out1;

        r[0] = false;
        circuit_evaluate(c, x, r, &out0);

        r[0] = true;
        circuit_evaluate(c, x, r, &out1);

        printf("  %d  %d |   %d  |   %d\n", x[0], x[1], out0, out1);
    }

    /* Topological analysis */
    size_t depth = circuit_compute_depth(c);
    bool is_dag = circuit_is_dag(c);
    printf("\nCircuit depth: %zu  DAG: %s\n\n", depth, is_dag ? "yes" : "no");

    circuit_destroy(c);

    /* Part 2: PARITY circuit construction */
    printf("--- Part 2: PARITY Circuit Family ---\n\n");

    for (size_t n = 2; n <= 8; n *= 2) {
        BooleanCircuit *parity = circuit_parity(n);
        if (parity) {
            size_t d = circuit_compute_depth(parity);
            printf("PARITY_%zu: %zu gates, depth %zu\n",
                   n, circuit_count_gates(parity), d);
            circuit_destroy(parity);
        }
    }

    /* Part 3: Shannon's lower bound */
    printf("\n--- Part 3: Shannon's Counting Lower Bound ---\n\n");

    for (size_t n = 2; n <= 10; n++) {
        size_t lb = shannon_lower_bound(n, 3, 2);
        printf("n=%2zu: Most functions need >= %6zu gates (out of %zu possible functions)\n",
               n, lb,
               (size_t)(pow(2.0, pow(2.0, (double)n)) / 1e6));
    }

    /* Part 4: Adleman's Theorem demonstration */
    printf("\n--- Part 4: Adleman's Theorem (BPP ⊆ P/poly) ---\n\n");

    for (size_t n = 1; n <= 20; n++) {
        double err = pow(2.0, -(double)(2 * n));  /* Error = 2^{-2n} */
        size_t advice_size = 0;
        bool exists = adleman_construction(n, err, &advice_size);

        printf("n=%2zu: 2^n=%7zu inputs, error=2^{-%zu}, ",
               n, (size_t)(1ULL << n), 2*n);
        if (exists) {
            printf("advice exists (%zu bits)\n", advice_size);
        } else {
            printf("union bound fails (>1)\n");
        }
    }

    /* Part 5: Concentration bounds comparison */
    printf("\n--- Part 5: Chernoff vs Hoeffding vs Exact ---\n\n");

    printf("Chernoff bound vs exact binomial tail:\n");
    printf("  n=100, p=0.5, delta=0.2:\n");
    double chernoff, exact;
    chernoff_vs_exact(100, 0.5, 0.2, &chernoff, &exact);
    printf("    Chernoff: %.6f\n", chernoff);
    printf("    Exact:    %.6f\n", exact);
    printf("    Ratio:    %.2f (Chernoff / Exact)\n", chernoff / exact);

    printf("\n  n=1000, p=0.5, delta=0.1:\n");
    chernoff_vs_exact(1000, 0.5, 0.1, &chernoff, &exact);
    printf("    Chernoff: %.6f\n", chernoff);
    printf("    Exact:    %.6f\n", exact);
    printf("    Ratio:    %.2f\n\n", chernoff / exact);

    printf("Hoeffding sample sizes for (ε,δ)-approximation:\n");
    double epsilons[] = {0.1, 0.05, 0.01};
    double deltas[] = {0.1, 0.05, 0.01};
    for (int ei = 0; ei < 3; ei++) {
        for (int di = 0; di < 3; di++) {
            size_t n = hoeffding_sample_size(epsilons[ei], deltas[di]);
            printf("  ε=%.2f δ=%.2f → n=%zu\n", epsilons[ei], deltas[di], n);
        }
    }

    printf("\n=== Demo Complete ===\n");
    return 0;
}
