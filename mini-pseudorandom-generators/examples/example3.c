#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "prg_core.h"
#include "prg_constructions.h"
#include "hardcore_predicates.h"
#include "next_bit_tests.h"
#include "nist_tests.h"

int main(void) {
    printf("=== PRG Example 3: Cryptographic Applications ===\n\n");

    printf("1. Hash_DRBG (NIST SP 800-90A):\n");
    uint8_t entropy[32], nonce[16];
    for (int i = 0; i < 32; i++) entropy[i] = (uint8_t)(i * 11 + 17);
    for (int i = 0; i < 16; i++) nonce[i] = (uint8_t)(i * 7 + 3);
    HashDRBG drbg = hash_drbg_instantiate(entropy, 32, nonce, 16, NULL, 0);
    uint8_t drbg_out[128];
    hash_drbg_generate(&drbg, drbg_out, 128, NULL, 0);
    printf("   Generated 128 bytes from Hash_DRBG\n");
    printf("   First 16 bytes: ");
    for (int i = 0; i < 16; i++) printf("%02x ", drbg_out[i]);
    printf("\n\n");

    printf("2. CTR_DRBG (AES-CTR based):\n");
    uint8_t seed[48];
    for (int i = 0; i < 48; i++) seed[i] = (uint8_t)(i + 1);
    CTRDRBG ctr = ctr_drbg_create(seed, 48);
    uint8_t ctr_out[64];
    ctr_drbg_generate(&ctr, ctr_out, 64);
    printf("   Generated 64 bytes from CTR_DRBG\n");
    printf("   First 16 bytes: ");
    for (int i = 0; i < 16; i++) printf("%02x ", ctr_out[i]);
    printf("\n\n");

    printf("3. ANSI X9.31 PRNG:\n");
    uint8_t x931_seed[16];
    for (int i = 0; i < 16; i++) x931_seed[i] = (uint8_t)(i * 255 / 16);
    uint8_t x931_out[64];
    ansi_x931_prng(x931_seed, 16, x931_out, 64, 20240601ULL);
    printf("   Generated 64 bytes from ANSI X9.31\n");
    printf("   First 16 bytes: ");
    for (int i = 0; i < 16; i++) printf("%02x ", x931_out[i]);
    printf("\n\n");

    printf("4. NIST Test Suite on PRG output:\n");
    bool test_bits[2048];
    for (int i = 0; i < 2048; i++)
        test_bits[i] = ((i * 25214903917ULL + 11) >> 16) & 1;

    NISTResult results[14];
    bool all_pass = nist_test_suite(test_bits, 2048, results, 14);

    const char *test_names[] = {
        "Monobit", "Block Freq", "Runs", "Longest Run",
        "Rank", "DFT", "Template", "Universal",
        "Linear Comp", "Serial", "Entropy", "Cum Sum",
        "Excursions", "Excursions Var"
    };
    for (int i = 0; i < 4; i++) {
        printf("   %-12s: P=%.6f [%s]\n",
               test_names[results[i].type],
               results[i].p_value,
               results[i].passed ? "PASS" : "FAIL");
    }
    printf("   All tests passed: %s\n\n", all_pass ? "yes" : "no");

    printf("5. Miller-Rabin Primality Test:\n");
    uint64_t primes[] = {2, 3, 7, 97, 104729};
    uint64_t composites[] = {1, 4, 100, 104730};
    for (int i = 0; i < 5; i++)
        printf("   %lu: %s\n", (unsigned long)primes[i],
               is_prime_u64(primes[i]) ? "prime" : "composite");
    for (int i = 0; i < 4; i++)
        printf("   %lu: %s\n", (unsigned long)composites[i],
               is_prime_u64(composites[i]) ? "prime" : "composite");
    printf("\n");

    printf("6. Modular Exponentiation:\n");
    printf("   2^100 mod 1000 = %lu\n", (unsigned long)modpow_u64(2, 100, 1000));
    printf("   3^50  mod 7    = %lu\n", (unsigned long)modpow_u64(3, 50, 7));
    printf("   5^20  mod 97   = %lu\n", (unsigned long)modpow_u64(5, 20, 97));
    printf("\n");

    printf("7. Blum-Blum-Shub Hardcore Bits:\n");
    bool bbs_bits[20];
    hc_blum_blum_shub(3, 11*19, 20, bbs_bits); /* N=209 */
    printf("   First 20 BBS bits: ");
    for (int i = 0; i < 20; i++) printf("%d", bbs_bits[i]);
    printf("\n\n");

    printf("8. Stream Encryption (PRG as stream cipher):\n");
    bool key[8] = {1,0,1,0,1,1,0,0};
    bool plaintext[16] = {0,1,0,1,1,0,0,1,1,1,0,0,1,0,1,1};
    bool ciphertext[16] = {0};
    PRG prg;
    prg_initialize(&prg, 8, 16, NULL, "DemoPRG");
    prg_stream_encrypt(key, 8, plaintext, 16, ciphertext, &prg);
    printf("   Plaintext:  ");
    for (int i = 0; i < 16; i++) printf("%d", plaintext[i]);
    printf("\n   With NULL gen, ct is unchanged.\n");

    printf("\n=== Example Complete ===\n");
    return 0;
}
