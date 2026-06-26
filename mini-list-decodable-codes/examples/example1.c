/*
 * example1.c — Reed-Solomon Code Encoding and List-Decoding Demo
 *
 * Demonstrates:
 *   1. Creating an RS code [7,3,5]_7 over GF(7)
 *   2. Encoding a message (polynomial evaluation)
 *   3. Corrupting the codeword with errors
 *   4. Welch-Berlekamp unique decoding
 *   5. Sudan list-decoding (beyond unique-decoding radius)
 *   6. Johnson bound verification
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "finite_field.h"
#include "rs_codes.h"
#include "list_decode_core.h"

int main(void) {
    printf("=== Reed-Solomon Code with List Decoding ===\n\n");

    /* Step 1: Create RS code [7,3,5]_7 over GF(7) */
    printf("Step 1: Creating RS code [7,3,5]_7 over GF(7)...\n");
    RSCode *rs = rs_create(7, 3, 7, NULL);
    printf("  n=%zu, k=%zu, d=%zu, q=%zu\n", rs->n, rs->k,
           rs_min_distance(rs), rs->q);

    /* Step 2: Encode message f(x) = 1 + 2x + 3x^2 */
    printf("\nStep 2: Encoding message f(x) = 1 + 2x + 3x^2...\n");
    size_t msg[] = {1, 2, 3};
    size_t codeword[7];
    rs_encode(rs, msg, 3, codeword);

    printf("  Codeword: [");
    for (size_t i = 0; i < 7; i++) printf("%zu ", codeword[i]);
    printf("]\n");

    printf("  Verifying: f(0)=%zu, f(1)=%zu\n", codeword[0], codeword[1]);

    /* Step 3: Corrupt with 1 error */
    printf("\nStep 3: Corrupting codeword with 1 error...\n");
    size_t received[7];
    memcpy(received, codeword, 7 * sizeof(size_t));
    received[2] = (received[2] + 3) % 7;
    printf("  Codeword[2] changed from %zu to %zu\n", codeword[2], received[2]);

    /* Step 4: Welch-Berlekamp unique decoder */
    printf("\nStep 4: Welch-Berlekamp unique decoder...\n");
    size_t decoded[3] = {0};
    bool wb_ok = rs_welch_berlekamp(rs, received, decoded);
    printf("  Decoding %s\n", wb_ok ? "succeeded" : "failed");
    if (wb_ok) {
        printf("  Decoded message: [%zu, %zu, %zu]\n",
               decoded[0], decoded[1], decoded[2]);
    }

    /* Step 5: Sudan list-decoding */
    printf("\nStep 5: Sudan list-decoding algorithm...\n");
    size_t *sudan_out = NULL;
    size_t L = rs_sudan_list_decode(rs, received, 1, &sudan_out);
    printf("  Found %zu candidate(s)\n", L);
    if (L > 0 && sudan_out) {
        printf("  Candidate 0 coefficients: [%zu, %zu, %zu]\n",
               sudan_out[0], sudan_out[1], sudan_out[2]);
        free(sudan_out);
    }

    /* Step 6: Johnson bound analysis */
    printf("\nStep 6: Johnson bound analysis for RS [7,3,5]_7:\n");
    double jr = rs_johnson_radius(rs);
    printf("  Johnson radius: %.2f (error count for poly-list-size)\n", jr);
    double cap = rs_list_capacity_limit(7, 3.0 / 7.0);
    printf("  List-decoding capacity: delta = %.4f\n", cap);

    printf("\n=== Example Complete ===\n");

    rs_free(rs);
    return 0;
}
