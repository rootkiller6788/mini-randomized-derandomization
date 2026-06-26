int main(void){
#include <stdio.h>
#include "nw_prg.h"
printf("Testing PRG construct optimal...\n");
NWPRG *prg = nw_prg_construct_optimal(50, 1000.0, NW_PRG_STANDARD); if(prg) printf("OK: k=%zu m=%zu\n", prg->seed_len, prg->output_len); else printf("NULL\n");
printf("Testing statistical tests...\n");
StatisticalTestResult *r = nw_prg_statistical_tests(prg, 5, 128); if(r) { printf("OK: tests=%zu passes=%d\n", r->num_tests, r->passes_all); nw_statistical_test_result_free(r); }
printf("Done.\n"); return 0; }
