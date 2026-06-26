/* demo_nw_prg.c — Interactive NW PRG Parameter Explorer */
#include <stdio.h>
#include <math.h>
#include "nw_core.h"
#include "nw_prg.h"
int main(void) {
    printf("NW PRG Parameter Explorer
");
    printf("Hardness S | Output m | Seed k | Stretch
");
    for (double S = 100; S <= 1000000; S *= 10) {
        for (size_t m = 16; m <= 1024; m *= 2) {
            size_t k = nw_seed_length(m, S);
            printf("  %.0e | %6zu | %6zu | %.1f
", S, m, k, (double)m/k);
        }
    }
    return 0;
}
