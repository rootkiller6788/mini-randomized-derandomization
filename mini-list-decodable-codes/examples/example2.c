/*
 * example2.c — Folded Reed-Solomon Codes: List-Decoding to Capacity
 *
 * Demonstrates:
 *   1. Creating a folded RS code (m=2, folding parameter)
 *   2. Encoding a message in the folded representation
 *   3. Computing the list-decoding radius for folded RS
 *   4. Verifying capacity-approaching property
 *   5. List-decoding via linear-algebraic approach
 *   6. Comparing with unfolded RS capacity
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "folded_rs.h"
#include "finite_field.h"

int main(void) {
    printf("=== Folded RS Codes: Achieving List-Decoding Capacity ===\n\n");

    /* Create folded RS code with n=8, k=4, q=11, folding m=2 */
    printf("Step 1: Creating folded RS code (n=8, k=4, q=11, m=2)...\n");
    FoldedRSCode *frs = frs_create(8, 4, 11, 2, NULL);
    if (!frs) {
        printf("  Failed to create folded RS code.\n");
        return 1;
    }
    printf("  Folded parameters: N=%zu, K=%zu, rate=%.2f\n",
           frs->N, frs->K, frs_rate(frs));

    /* Encode a message */
    printf("\nStep 2: Encoding message [1, 2, 3, 4]...\n");
    size_t msg[] = {1, 2, 3, 4};
    size_t *codeword = NULL;
    frs_encode(frs, msg, &codeword);
    if (codeword) {
        printf("  Codeword (unfolded, %zu symbols): [", frs->n);
        for (size_t i = 0; i < frs->n; i++) printf("%zu ", codeword[i]);
        printf("]\n");
        printf("  Folded representation: %zu blocks of %zu symbols each\n",
               frs->N, frs->m);
        free(codeword);
    }

    /* List-decoding radius analysis */
    printf("\nStep 3: List-decoding radius analysis...\n");
    double radius = frs_list_decoding_radius(frs);
    double rate = frs_rate(frs);
    double capacity = 1.0 - rate;
    printf("  Rate R = %.2f, Capacity delta_cap = %.4f\n", rate, capacity);
    printf("  Achievable radius (m=%zu): %.4f\n", frs->m, radius);
    printf("  Gap to capacity: %.4f\n", capacity - radius);

    /* Capacity-approaching verification */
    printf("\nStep 4: Capacity-approaching verification...\n");
    bool cap_ok = frs_verify_capacity_approaching(frs, 0.3);
    printf("  Within epsilon=0.30 of capacity: %s\n",
           cap_ok ? "YES" : "NO");

    /* Effect of larger folding parameter */
    printf("\nStep 5: Larger folding parameters...\n");
    for (size_t m_val = 2; m_val <= 8; m_val *= 2) {
        FoldedRSCode *frs_m = frs_create(8 * m_val, 4 * m_val, 11, m_val, NULL);
        if (frs_m) {
            double r_m = frs_list_decoding_radius(frs_m);
            double gap = (1.0 - frs_rate(frs_m)) - r_m;
            printf("  m=%zu: radius=%.4f, gap=%.4f\n", m_val, r_m, gap);
            frs_free(frs_m);
        }
    }

    /* Theorem: Guruswami-Rudra (2008) */
    printf("\nStep 6: Guruswami-Rudra Theorem (2008):\n");
    printf("  \"For any epsilon > 0, folded RS codes with\n");
    printf("   sufficiently large folding parameter m\n");
    printf("   achieve list-decoding radius delta = 1-R-epsilon.\"\n");
    printf("  As m -> infinity, the gap approaches 0.\n");

    printf("\n=== Example Complete ===\n");

    frs_free(frs);
    return 0;
}
