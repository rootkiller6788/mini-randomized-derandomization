/************************************************************
 * example3.c -- Zig-Zag Product & SL Connectivity Demo
 *
 * Demonstrates the zig-zag product construction (Reingold-Vadhan-
 * Wigderson 2002) and the SL = L connectivity algorithm.
 *
 * The zig-zag product transforms a large-degree expander G and a
 * small expander H into a new expander G(z)H with dramatically
 * reduced degree while preserving expansion. This is the key
 * building block in Reingold's proof that SL = L.
 *
 * Build: make example3
 ************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "expander_core.h"
#include "expander_constructions.h"
#include "expander_applications.h"
#include "zigzag_product.h"

int main(void) {
    printf("=== Example 3: Zig-Zag Product & SL Connectivity ===\n\n");

    /* Step 1: Build base graphs */
    printf("Step 1: Building base graphs...\n");

    /* G: 8-regular Margulis expander on 9 vertices (3x3 grid) */
    ExpanderGraph *G = exp_margulis_gabber_galil(9);
    if (!G) { printf("G construction failed\n"); return 1; }
    printf("  G: n=%zu, d=%zu (Margulis)\n", G->n, G->d);

    /* H: Must have n_H = d_G = 8 vertices.
     * Build a small expander on 8 vertices using random construction */
    ExpanderGraph *H = exp_random_regular(8, 3);
    if (!H) { printf("H construction failed\n"); exp_free(G); return 1; }
    printf("  H: n=%zu, d=%zu (Random regular)\n", H->n, H->d);

    /* Step 2: Compute zig-zag product */
    printf("\nStep 2: Computing zig-zag product G(z)H...\n");
    ZigZagResult zr = zz_compute(G, H);
    if (!zr.result) {
        printf("Zig-zag product failed (n_H must equal d_G = %zu)\n", G->d);
        exp_free(G); exp_free(H);
        return 1;
    }
    printf("  Result: n=%zu, d=%zu\n", zr.result->n, zr.result->d);

    /* Verify construction */
    if (zz_verify_construction(&zr))
        printf("  Construction integrity: VERIFIED\n");

    /* Step 3: Spectral analysis */
    printf("\nStep 3: Spectral Analysis...\n");
    double lamG = exp_second_eigenvalue(G);
    double lamH = exp_second_eigenvalue(H);
    double lamZ = exp_second_eigenvalue(zr.result);

    printf("  lambda(G):      %.6f\n", lamG);
    printf("  lambda(H):      %.6f\n", lamH);
    printf("  lambda(G(z)H):  %.6f\n", lamZ);
    printf("  Theoretical bound: %.6f\n", zz_spectral_bound(G, H));
    printf("  Empirical/theoretical ratio: %.4f\n",
           zz_empirical_lambda(&zr));

    /* Degree reduction */
    double degree_ratio = zz_degree_reduction(G->d, H->d);
    printf("  Degree reduction: %.4fx (d_H^2/d_G = %zu/%zu)\n",
           degree_ratio, H->d * H->d, G->d);

    /* Step 4: SL Connectivity demonstration */
    printf("\nStep 4: SL Connectivity (Reingold's algorithm)...\n");

    /* Test connectivity on the expander */
    bool connected;
    exp_sl_connectivity(zr.result, 0, zr.result->n / 2, &connected);
    printf("  s=0 to t=%zu: %s\n", zr.result->n / 2,
           connected ? "CONNECTED" : "NOT CONNECTED");

    /* Path length estimate */
    size_t path_len = exp_sl_path_length(zr.result, 0, zr.result->n / 2);
    printf("  Estimated path length: %zu\n", path_len);

    /* Step 5: Derandomized sampling on the zig-zag product */
    printf("\nStep 5: Derandomized Sampling...\n");
    size_t samples[50];
    exp_derandomized_sampling(zr.result, 50, samples);
    double quality = exp_sampling_quality(samples, 50, zr.result->n);
    printf("  Sampling quality (TV distance): %.6f\n", quality);
    printf("  (0 = perfect uniform, 1 = worst case)\n");

    /* Step 6: Expander code demonstration */
    printf("\nStep 6: Expander Code Demo...\n");
    char message[] = "11010010";
    char codeword[32];
    printf("  Original message: %s\n", message);
    exp_code_encode(message, 8, zr.result, codeword);
    printf("  Encoded (msg+parity): %s\n", codeword);

    /* Corrupt one bit and decode */
    codeword[2] = (codeword[2] == '0') ? '1' : '0';
    printf("  Corrupted bit 2:    %s\n", codeword);
    exp_code_flip_decode(codeword, 8, zr.result, 5);
    printf("  After decoding:     %s\n", codeword);

    /* Cleanup */
    zz_free_result(&zr);
    exp_free(G);
    exp_free(H);

    printf("\n=== Example 3 Complete ===\n");
    return 0;
}
