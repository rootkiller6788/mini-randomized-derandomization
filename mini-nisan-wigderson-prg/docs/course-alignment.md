# Course Alignment — mini-nisan-wigderson-prg

## Nine-School Curriculum Mapping

| School | Course | Relevant Chapters | Topics Covered |
|--------|--------|-------------------|----------------|
| **MIT** | 6.841 Advanced Complexity | Ch 20: Derandomization | NW PRG construction, hardness amplification, IW97 |
| **Stanford** | CS254 Computational Complexity | Ch 19-20: Pseudorandomness | PRG definitions, hybrid argument, derandomization |
| **Berkeley** | CS278 Computational Complexity | Ch 11: Derandomization | NW Theorem, BPP derandomization, expanders |
| **CMU** | 15-855 Graduate Complexity | Lectures 18-22 | Design constructions, hardness vs randomness |
| **Princeton** | COS 522 Computational Complexity | Ch 8: Randomness | Cryptographic PRGs, hard-core bits, GL theorem |
| **Caltech** | CS 151 Complexity Theory | Ch 14: Derandomization | Circuit lower bounds, switching lemma, AC0 |
| **Cambridge** | Part II Complexity Theory | Ch 7: Randomised Computation | BPP, RP, derandomization |
| **Oxford** | Advanced Complexity Theory | Ch 6: Pseudorandomness | NW generator, hardness assumptions |
| **ETH** | 263-4650 Advanced Complexity | Ch 8-9: Circuit Complexity | AC0 lower bounds, Hastad, Razborov-Smolensky |

## Specific Course Module Alignment

### MIT 6.841 — Advanced Complexity Theory

| Module | Our Implementation |
|--------|-------------------|
| 20.1: Pseudorandom generators | `nw_prg.h/c` — NW PRG construction/evaluation |
| 20.2: Hardness vs randomness | `nw_core.h` — Core definitions |
| 20.2.1: Combinatorial designs | `nw_designs.h/c` — RS, random, expander designs |
| 20.2.2: NW construction | `nw_prg.c` — `nw_prg_construct_optimal` |
| 20.3: Derandomizing BPP | `nw_prg.c` — `nw_bpp_derandomize` |
| 20.4: IW97 (BPP=P) | `nw_hardness.c` — `nw_iw_derandomize` |

### Stanford CS254 — Computational Complexity

| Module | Our Implementation |
|--------|-------------------|
| 19 Yao XOR Lemma | `nw_hardness.c` — `nw_yao_xor_amplify` |
| 19 Goldreich-Levin | `nw_circuits.c` — `nw_goldreich_levin_hard_predicate` |
| 20 NW PRG | `nw_prg.c` — full pipeline |
| 20 Hybrid argument | `nw_prg.c` — `nw_hybrid_position` |

### Berkeley CS278 — Computational Complexity

| Module | Our Implementation |
|--------|-------------------|
| 11 Circuit complexity | `nw_circuits.c` — Shannon, Lupanov bounds |
| 11 Switching lemma | `nw_circuits.c` — `nw_hastad_applies` |
| 11 AC0 lower bounds | `nw_circuits.c` — `nw_ac0_parity_lower_bound` |
| 11 Natural proofs | `nw_hardness.c` — `nw_natural_proofs_barrier` |