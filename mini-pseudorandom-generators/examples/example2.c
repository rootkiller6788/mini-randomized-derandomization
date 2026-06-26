#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "prg_core.h"
#include "prg_constructions.h"
#include "hardcore_predicates.h"

int main(void) {
    printf("=== PRG Example 2: Construction and Statistical Testing ===\n\n");

    printf("1. LCG Implementation:\n");
    LCG lcg = lcg_create(42, 1103515245, 12345, 2147483648);
    printf("   First 5 values from LCG:\n");
    for (int i = 0; i < 5; i++) {
        printf("     X_%d = %lu\n", i, (unsigned long)lcg_next(&lcg));
    }
    double q = lcg_quality_test(&lcg, 10000);
    printf("   Quality score: %.4f\n\n", q);

    printf("2. Mersenne Twister MT19937:\n");
    MT19937 mt = mt19937_create(5489ULL);
    printf("   First 5 values from MT19937:\n");
    for (int i = 0; i < 5; i++) {
        printf("     Y_%d = %lu\n", i, (unsigned long)mt19937_extract(&mt));
    }
    bool bits[1000];
    mt = mt19937_create(12345ULL);
    mt19937_generate_bits(&mt, bits, 1000);
    double h = binary_entropy(bits, 1000);
    printf("   Binary entropy of 1000 bits: %.6f (ideal: 1.0)\n\n", h);

    printf("3. Xorshift128:\n");
    Xorshift128 xs = xorshift128_create(12345, 67890);
    printf("   First 5 values from Xorshift128:\n");
    for (int i = 0; i < 5; i++) {
        printf("     Z_%d = %lu\n", i, (unsigned long)xorshift128_next(&xs));
    }
    printf("\n");

    printf("4. ChaCha20 Stream Cipher:\n");
    uint8_t key[32], nonce[12];
    for (int i = 0; i < 32; i++) key[i] = (uint8_t)(i * 3 + 7);
    for (int i = 0; i < 12; i++) nonce[i] = (uint8_t)(i * 5 + 13);
    ChaCha20PRG chacha = chacha20_init(key, nonce, 0);
    uint8_t output[64];
    chacha20_keystream(&chacha, output, 64);
    printf("   First 16 bytes of keystream: ");
    for (int i = 0; i < 16; i++) printf("%02x ", output[i]);
    printf("\n\n");

    printf("5. Statistical Test Battery:\n");
    bool test_bits[1024];
    for (int i = 0; i < 1024; i++)
        test_bits[i] = ((i * 6364136223846793005ULL) >> 32) & 1;

    double p_mono = nist_frequency_monobit(test_bits, 1024);
    printf("   Monobit P-value:       %.6f\n", p_mono);

    double p_runs = nist_runs_test(test_bits, 1024);
    printf("   Runs P-value:          %.6f\n", p_runs);

    if (1024 >= 128) {
        double p_lr = nist_longest_run_ones(test_bits, 1024);
        printf("   Longest Run P-value:   %.6f\n", p_lr);
    }

    double p_cs = nist_cumulative_sums(test_bits, 1024, true);
    printf("   Cumulative Sum P-value: %.6f\n", p_cs);

    size_t lc = linear_complexity(test_bits, 1024);
    printf("   Linear complexity:     %zu (expected ~ %zu)\n", lc, (size_t)512);

    double ac = autocorrelation(test_bits, 1024, 1);
    printf("   Autocorrelation(lag=1): %.6f\n", ac);
    printf("\n");

    printf("6. Hardcore Predicate: Goldreich-Levin\n");
    bool x[64], r[64], result;
    for (int i = 0; i < 64; i++) {
        x[i] = (i * 7 + 3) & 1;
        r[i] = (i * 13 + 5) & 1;
    }
    hc_goldreich_levin(x, r, 64, &result);
    printf("   <x,r> mod 2 = %d\n", result);

    printf("\n=== Example Complete ===\n");
    return 0;
}
