/**
 * example3.c — Circuit Lower Bounds and Hardness
 *
 * Demonstrates:
 *   L4: Shannon's counting lower bound
 *   L6: PARITY and MAJORITY circuit constructions
 *   L8: Hastad switching lemma analysis
 *   L5: Goldreich-Levin hard predicate
 *
 * This example explores the circuit complexity side
 * of the hardness vs randomness paradigm.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nw_core.h"
#include "nw_circuits.h"
#include "nw_hardness.h"

int main(void) {
    printf("╔══════════════════════════════════════════╗\n");
    printf("║  Circuit Complexity & Hardness Analysis ║\n");
    printf("╚══════════════════════════════════════════╝\n\n");

    /* ── Shannon's Lower Bound ── */
    printf("── Shannon's Counting Lower Bound ────────\n");
    printf("n | Lower Bound (gates) | Upper Bound\n");
    printf("──┼──────────────────────┼────────────\n");
    for (size_t n = 2; n <= 16; n *= 2) {
        double lb = nw_shannon_lower_bound(n);
        double ub = nw_lupanov_upper_bound(n);
        printf("%2zu | %20.1f | %11.1f\n", n, lb, ub);
    }
    printf("\n");

    /* ── PARITY Circuit ── */
    printf("── PARITY Circuit Construction ───────────\n");
    for (size_t n = 2; n <= 8; n *= 2) {
        Circuit *c = nw_circuit_parity(n);
        if (c) {
            size_t gates = nw_circuit_gate_count(c);
            size_t depth = nw_circuit_depth(c);
            printf("PARITY-%zu: gates=%zu, depth=%zu\n", n, gates, depth);

            /* Verify on all inputs */
            bool input[16] = {0};
            size_t total = (size_t)1 << n;
            for (size_t i = 0; i < total; i++) {
                for (size_t j = 0; j < n; j++)
                    input[j] = (i >> j) & 1;
                bool expected = false;
                for (size_t j = 0; j < n; j++)
                    expected ^= input[j];
                bool actual = nw_circuit_eval_single(c, input);
                assert(expected == actual);
            }
            printf("  All %zu inputs verified ✓\n", total);
            nw_circuit_free(c);
        }
    }
    printf("\n");

    /* ── MAJORITY Circuit ── */
    printf("── MAJORITY Circuit Construction ─────────\n");
    for (size_t n = 3; n <= 7; n += 2) {
        Circuit *c = nw_circuit_majority(n);
        if (c) {
            printf("MAJORITY-%zu: gates=%zu\n", n, nw_circuit_gate_count(c));
            nw_circuit_free(c);
        }
    }
    printf("\n");

    /* ── AC0 Lower Bounds via Switching Lemma ── */
    printf("── AC0 Lower Bounds (Hastad) ─────────────\n");
    printf("d | n=100   | n=1000  | n=10000\n");
    printf("──┼──────────┼─────────┼────────\n");
    for (size_t d = 2; d <= 5; d++) {
        printf("%zu | ", d);
        for (size_t n = 100; n <= 10000; n *= 10) {
            bool applies = nw_hastad_applies(d, n);
            double lb = nw_ac0_parity_lower_bound(n, d);
            if (applies)
                printf("%8.1f | ", lb);
            else
                printf("  N/A    | ");
        }
        printf("\n");
    }
    printf("\n");

    /* ── Goldreich-Levin Hard Predicate ── */
    printf("── Goldreich-Levin Hard Predicate ────────\n");
    Circuit *gl = nw_goldreich_levin_hard_predicate(4);
    if (gl) {
        printf("GL-4: inputs=%zu, gates=%zu\n",
               gl->n_inputs, nw_circuit_gate_count(gl));

        /* Test: IP(x,r) = XOR(x_i AND r_i) */
        bool test_in[8] = {true, false, true, false,  /* x */
                           false, true, true, false}; /* r */
        /* x·r = 1·0 + 0·1 + 1·1 + 0·0 = 1 → odd parity → true */
        bool result = nw_circuit_eval_single(gl, test_in);
        printf("IP((1010),(0110)) = %s (expected: true)\n",
               result ? "true" : "false");

        nw_circuit_free(gl);
    }
    printf("\n");

    /* ── Natural Proofs Barrier ── */
    printf("── Natural Proofs Barrier ────────────────\n");
    for (size_t n = 8; n <= 32; n *= 2) {
        double thresh = nw_razborov_rudich_threshold(n);
        printf("n=%2zu: Razborov-Rudich threshold ≈ 2^%.0f, ",
               n, log2(thresh));
        bool barrier = nw_natural_proofs_barrier(thresh * 2, n);
        printf("barrier=%s\n", barrier ? "BREACHED" : "safe");
    }

    printf("\n── Done ──────────────────────────────────\n");
    return 0;
}
