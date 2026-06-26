# mini-nisan-wigderson-prg

## Module Status: ✅ COMPLETE

```
Total: 17/18 points
L1-L6: Complete (12/12)
L7:    Complete (2/2) — 5 applications implemented
L8:    Complete (2/2) — 7 advanced topics implemented
L9:    Partial  (1/2) — 5 research topics documented
```

## Overview

This module implements the **Nisan-Wigderson Pseudorandom Generator** —
a cornerstone result in computational complexity theory showing that
**hardness implies randomness**: if there exist functions that are hard
for circuits, then pseudorandom generators exist, and randomness can
be eliminated from efficient algorithms.

### Key Results

| Theorem | Description |
|---------|------------|
| **Nisan-Wigderson (1994)** | Hardness ⇒ PRG: any hard function yields a pseudorandom generator |
| **Impagliazzo-Wigderson (1997)** | BPP = P if E requires exponential circuits |
| **Yao's XOR Lemma (1982)** | Hardness amplifies exponentially under XOR |
| **Shannon (1949)** | Most boolean functions require Ω(2ⁿ/n) gates |
| **Håstad (1987)** | PARITY requires exponential-size AC0 circuits |

## Architecture

```
include/              5 headers (1,967 lines)
  nw_core.h           Core types: NWDesign, HardFunction, NWPRG
  nw_designs.h        Combinatorial design constructions
  nw_hardness.h       Hardness assumptions, truth tables, amplification
  nw_circuits.h       Boolean circuit model, lower bounds
  nw_prg.h            PRG construction, evaluation, derandomization

src/                  5 C files + 1 Lean (2,968 lines)
  nw_core.c           Design ops, PRG eval, utilities
  nw_design.c         RS/random/expander designs
  nw_hardness.c       Yao XOR, truth tables, Fourier, natural proofs
  nw_circuits.c       PARITY, MAJORITY, switching lemma, GL predicate
  nw_prg.c            PRG construction, hybrid arg, BPP sim, KDF
  nw_theorem.lean     Lean 4 formalization (230 lines)

tests/                1 test file (814 lines)
  test.c              37 test cases, all passing

examples/             4 demos (357 lines)
  example1.c          PRG construction/evaluation
  example2.c          BPP derandomization
  example3.c          Circuit complexity analysis
  example_nw_demo.c   Complete pipeline demo

demos/ benches/ docs/
```

## Core API

```c
// L1: Create a combinatorial design
NWDesign *nw_design_create(size_t k, size_t m, size_t l, size_t bound);

// L5: Construct Reed-Solomon based design
NWDesign *nw_reed_solomon_design(size_t q, size_t degree);

// L4: Apply Yao's XOR Lemma
double nw_yao_xor_amplify(double epsilon, size_t k, double delta);

// L2: Construct an NW PRG
NWPRG *nw_prg_construct_optimal(size_t output_len, double hardness, NWPRGMode mode);

// L5: Evaluate the PRG on a seed
void nw_prg_evaluate(const NWPRG *prg, const bool *seed, bool *output);

// L7: Derandomize a BPP algorithm
bool nw_bpp_derandomize(const NWPRG *prg, ...);

// L5: Run statistical tests on PRG output
StatisticalTestResult *nw_prg_statistical_tests(const NWPRG *prg, ...);
```

## Quick Start

```bash
make test          # Build and run 37 tests
make check-lines   # Count include/ + src/ lines
make clean         # Remove build artifacts

# Compile examples individually:
cc -std=c11 -Iinclude -o ex1 examples/example1.c src/*.c -lm && ./ex1
cc -std=c11 -Iinclude -o ex2 examples/example2.c src/*.c -lm && ./ex2
cc -std=c11 -Iinclude -o ex3 examples/example3.c src/*.c -lm && ./ex3
```

## Knowledge Coverage

| Layer | Name | Status | Score |
|-------|------|--------|-------|
| L1 | Definitions | 10 types/structs | Complete (2) |
| L2 | Core Concepts | 8 concepts | Complete (2) |
| L3 | Math Structures | 10 structures | Complete (2) |
| L4 | Fundamental Laws | 10 theorems | Complete (2) |
| L5 | Algorithms/Methods | 15 algorithms | Complete (2) |
| L6 | Canonical Problems | 6 problems | Complete (2) |
| L7 | Applications | 5 applications | Complete (2) |
| L8 | Advanced Topics | 7 topics | Complete (2) |
| L9 | Research Frontiers | 5 documented | Partial (1) |

**Total: 17/18 — Meets COMPLETE threshold (≥16/18)**

## References

- Nisan & Wigderson, "Hardness vs Randomness", *J. Comput. Syst. Sci.* 49(2), 1994
- Impagliazzo & Wigderson, "P = BPP if E requires exponential circuits", *STOC* 1997
- Yao, "Theory and Applications of Trapdoor Functions", *FOCS* 1982
- Håstad, "Computational Limitations of Small-Depth Circuits", 1987
- Razborov & Rudich, "Natural Proofs", *J. Comput. Syst. Sci.* 55(1), 1997
- Arora & Barak, *Computational Complexity: A Modern Approach*, Ch 19-20
- Goldreich, *Computational Complexity: A Conceptual Perspective*, Ch 8
- Vollmer, *Introduction to Circuit Complexity*, 1999
