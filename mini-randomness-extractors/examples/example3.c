/**
 * example3.c -- Applications of Randomness Extractors
 *
 * Demonstrates L7:
 *   1. BPP derandomization via seed enumeration
 *   2. Cryptographic key derivation from weak secrets
 *   3. Privacy amplification protocol
 *   4. Network extractor simulation
 *   5. Security bound computation
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "extractor_core.h"
#include "extractor_constructions.h"
#include "extractor_applications.h"

/* A simple BPP algorithm for demonstration:
 *   returns true with probability biased by seed bits
 */
static bool sample_bpp(const char *inp, size_t ilen,
                       const bool *seed, size_t slen) {
    (void)inp; (void)ilen;
    /* Majority of seed bits */
    size_t ones = 0;
    for (size_t i = 0; i < slen; i++) if (seed[i]) ones++;
    return ones > slen / 2;
}

int main(void) {
    printf("============================================================\n");
    printf("  Example 3: Applications\n");
    printf("============================================================\n\n");

    /* Part 1: BPP derandomization */
    printf("--- Part 1: BPP Derandomization ---\n");
    {
        bool result = false;
        bool ok = ext_derandomize_bpp(sample_bpp, "input", 5, 8, &result);
        printf("  Algorithm: majority-of-seed-bits\n");
        printf("  Enumerated 2^8 = 256 seeds\n");
        printf("  Derandomized result: %s\n", result ? "ACCEPT" : "REJECT");
        printf("  Majority decisive: %s\n\n", ok ? "YES" : "NO");
    }

    /* Part 2: Cryptographic key derivation */
    printf("--- Part 2: Key Derivation ---\n");
    {
        uint8_t dh_secret[8] = {
            0xAB, 0xCD, 0xEF, 0x01, 0x23, 0x45, 0x67, 0x89
        };
        CryptoKey key = ext_key_derivation(dh_secret, 8, 32);
        printf("  DH shared secret: ABCDEF0123456789\n");
        printf("  Derived key:       ");
        for (size_t i = 0; i < 8 && i < key.len; i++)
            printf("%02X", key.key[i]);
        printf("... (%zu bytes)\n", key.len);

        double eps = ext_key_derivation_bound(64, 48.0, 256);
        printf("  Security: eps <= 2^{-%d} = %.2e\n",
               (int)((48.0-32.0)/2.0), eps);
        ck_free(&key);
        printf("\n");
    }

    /* Part 3: Privacy amplification */
    printf("--- Part 3: Privacy Amplification ---\n");
    {
        uint8_t shared[8] = {0x01,0x23,0x45,0x67,0x89,0xAB,0xCD,0xEF};
        uint8_t seed[4]   = {0xCA,0xFE,0xBA,0xBE};
        uint8_t key[8]    = {0};

        bool ok = ext_privacy_amplification(shared, 8, seed, 4, key, 8);
        printf("  Shared secret (weak): 0123456789ABCDEF\n");
        printf("  Public seed:          CAFEBABE\n");
        printf("  Amplified key:        ");
        for (int i = 0; i < 8; i++) printf("%02X", key[i]);
        printf("\n");
        printf("  Protocol success: %s\n", ok ? "YES" : "NO");

        double eps = ext_privacy_amplification_bound(64.0, 64);
        printf("  Security bound: eps = 2^{-%d/2}\n",
               (int)(64.0 - 64.0));
        printf("  (LHL: key is epsilon-close to uniform given Eve's view)\n\n");
    }

    /* Part 4: Network extractor */
    printf("--- Part 4: Network Extractor ---\n");
    {
        printf("  Scenario: 3 parties, each with 4-bit weak source\n");
        printf("  Adversary controls at most 1 party\n");
        printf("  Goal: extract 2-bit common random string\n\n");

        bool ok = ext_network_extractor_simulate(4, 3, 4, 2, 50);
        printf("  Network protocol success rate: %s\n",
               ok ? ">= 95%" : "< 95%");

        bool src1[4] = {1,0,1,0}, src2[4] = {0,1,1,0};
        bool src3[4] = {1,1,0,0}, src4[4] = {0,0,1,1};
        const bool *sources[4] = {src1, src2, src3, src4};
        bool output[2] = {0};
        /* Fault tolerant: even if 1 of 4 is corrupted */
        bool ft_ok = ext_multi_source_fault_tolerant(
            sources, 4, 4, 1, 4, 2, output);
        printf("  Fault-tolerant extraction: %s, output=%d%d\n\n",
               ft_ok ? "OK" : "FAIL", output[1], output[0]);
    }

    /* Part 5: Security bounds summary */
    printf("--- Part 5: Security Bounds Summary ---\n");
    {
        printf("  Classical LHL:       eps <= 2^{-(k-m)/2}\n");
        double q_eps = ext_quantum_proof_bound(128.0, 64);
        printf("  Quantum-proof:       eps <= 2 * 2^{-(k-m)/4}\n");
        printf("                       For k=128, m=64: eps_q = %.2e\n", q_eps);
        printf("\n  Theorem (Tomamichel et al., 2011):\n");
        printf("  Universal hashing is quantum-proof with\n");
        printf("  only a quadratic loss in parameters.\n");
    }

    printf("\n=== Example 3 Complete ===\n");
    return 0;
}
