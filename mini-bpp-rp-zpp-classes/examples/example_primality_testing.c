/******************************************************************************
 * example_primality_testing.c — ZPP/RP Primality Testing Demo
 *
 * Demonstrates Miller-Rabin primality testing as a ZPP (zero-error
 * probabilistic polynomial time) algorithm.
 *
 * Algorithm (Miller-Rabin):
 *   1. Write n-1 = 2^s * d (d odd)
 *   2. Pick random a ∈ [2, n-2]
 *   3. Compute x = a^d mod n
 *   4. If x = 1 or x = n-1: probably prime
 *   5. Square s-1 times: if any = n-1: probably prime
 *   6. Otherwise: composite (witness found)
 *
 * Reference: Miller (1976), Rabin (1980)
 * L5: Algorithms — randomized primality testing
 ******************************************************************************/

#include "randomized_classes.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

int main(void) {
    printf("=== ZPP Primality Testing: Miller-Rabin Demo ===\n\n");

    uint64_t test_numbers[] = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29,  /* small primes */
        4, 6, 8, 9, 10, 12, 14, 15, 16, 18,  /* small composites */
        97, 101, 103, 107, 109, 113,           /* medium primes */
        100, 121, 143, 169, 187,               /* composites: 10², 11², 11·13, 13², 11·17 */
        9973, 10007, 99991,                    /* larger primes */
        65521,                                  /* 2^16+1? No: 65521 is prime */
    };

    int num_tests = sizeof(test_numbers) / sizeof(test_numbers[0]);

    printf("Testing %d numbers with Miller-Rabin (randomized mode):\n\n", num_tests);
    printf("%-10s | %-10s | %-10s | %s\n", "Number", "Random", "Det", "Factors");
    printf("-----------+------------+------------+-------------------\n");

    RandomSource rs = random_source_init((uint64_t)time(NULL));

    int correct = 0;

    for (int i = 0; i < num_tests; i++) {
        uint64_t n = test_numbers[i];
        int result_rand = zpp_primality_test(n, &rs, false);
        int result_det = zpp_primality_test(n, &rs, true);

        /* Determine true primality by trial division */
        bool is_prime = true;
        if (n < 2) is_prime = false;
        else if (n == 2) is_prime = true;
        else if (n % 2 == 0) is_prime = false;
        else {
            for (uint64_t d = 3; d * d <= n && d < 10000; d += 2) {
                if (n % d == 0) { is_prime = false; break; }
            }
        }

        bool rand_correct = ((result_rand == 1) == is_prime);
        bool det_correct = ((result_det == 1) == is_prime);

        if (rand_correct) correct++;

        printf("%-10llu | %-10s | %-10s | ",
               (unsigned long long)n,
               rand_correct ? "OK" : "WRONG",
               det_correct ? "OK" : "WRONG");

        if (!is_prime) {
            /* Print small factors */
            printf("composite");
            for (uint64_t d = 2; d * d <= n && d < 50; d++) {
                if (n % d == 0) {
                    printf(" (%llu)", (unsigned long long)d);
                    if (n / d <= 50) printf("*%llu", (unsigned long long)(n/d));
                }
            }
        } else {
            printf("prime");
        }
        printf("\n");
    }

    printf("\nRandomized mode: %d/%d correct\n", correct, num_tests);

    /* Performance comparison */
    printf("\nPerformance test: testing numbers around 10^6:\n");
    uint64_t perf_nums[] = {1000003, 1000033, 1000037, 1000039, 2000003};
    for (int i = 0; i < 5; i++) {
        RandomSource rs2 = random_source_init(i * 31337ULL);
        int r = zpp_primality_test(perf_nums[i], &rs2, false);
        printf("  %llu: %s\n",
               (unsigned long long)perf_nums[i],
               r == 1 ? "prime" : (r == 0 ? "composite" : "unknown"));
    }

    printf("\nError analysis (Miller-Rabin):\n");
    printf("  Each random base: error <= 1/4 (probability of declaring\n");
    printf("  a composite as prime). With k random bases: error <= (1/4)^k.\n");
    printf("\n  k=10:  error <= 2^{-20} ~= 10^{-6}\n");
    printf("  k=20:  error <= 2^{-40} ~= 10^{-12}\n");
    printf("  k=40:  error <= 2^{-80} ~= 10^{-24}\n");

    return 0;
}
