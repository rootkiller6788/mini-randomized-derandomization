# mini-pseudorandom-generators

## Module Status: COMPLETE

**include/ + src/ = 3602 lines** (>= 3000) · **54/54 tests passing**

### Knowledge Coverage

| Level | Status | Score | Details |
|-------|--------|-------|---------|
| **L1** Definitions | **Complete** | 2 | PRG, Distinguisher, NextBitPredictor, LCG, MT19937, Xorshift128, ChaCha20, CTR_DRBG, Hash_DRBG, OneWayFn, HardcorePred |
| **L2** Core Concepts | **Complete** | 2 | Statistical distance, distinguish advantage, stretch/ratio, negligibility, computational distance, security level |
| **L3** Math Structures | **Complete** | 2 | GF(2) arithmetic, modular arithmetic, LCG period analysis, polynomial stretch composition, NW generator design, hybrid argument |
| **L4** Fundamental Theorems | **Complete** | 2 | Yao Theorem (1982), Hybrid Argument, Goldreich-Levin (1989), Blum-Micali (1984), BBS (1986), IW97, quantitative Yao |
| **L5** Algorithms/Methods | **Complete** | 2 | NIST SP 800-22 tests, cycle detection (Floyd), Berlekamp-Massey, BPP error reduction, stream encryption, NW generator, Miller-Rabin |
| **L6** Canonical Problems | **Complete** | 2 | PRG distinguisher problem, next-bit prediction, hardcore predicate extraction, one-way function inversion |
| **L7** Applications | **Complete** | 2 | Stream cipher (IND-CPA), BPP derandomization, NIST SP 800-90A DRBG (Hash_DRBG, CTR_DRBG), ANSI X9.31 |
| **L8** Advanced Topics | **Partial** | 1 | ChaCha20 (RFC 8439), Nisan-Wigderson framework, IW97 conditional derandomization, GL list decoding |
| **L9** Research Frontiers | **Partial** | 1 | Levin universal OWF, Impagliazzo hardcore lemma, Yao XOR lemma |

**Total: 16/18 -- COMPLETE**

### Core Definitions

- **PRG**: G:{0,1}^n -> {0,1}^{l(n)}, l(n)>n, PPT D has negligible advantage
- **Distinguisher Advantage**: adv(D,G,n) = |Pr[D(G(U_n))=1] - Pr[D(U_l)=1]|
- **Stretch**: s(n) = l(n) - n
- **Next-bit Unpredictability**: For all i, Pr[A(G(x)_{1..i}) = G(x)_{i+1}] <= 1/2 + negl(n)
- **Hardcore Predicate**: b(x) that is hard to predict given f(x) for OWF f

### Core Theorems

1. **Yao Theorem (1982)**: Next-bit unpredictability iff pseudorandomness
   - eps_dist <= m * eps_nbp (hybrid argument)
2. **Goldreich-Levin (1989)**: b(x,r)=<x,r> mod 2 is a universal hardcore predicate
3. **Blum-Micali (1984)**: MSB is hardcore for discrete log
4. **Blum-Blum-Shub (1986)**: LSB is hardcore for quadratic residues
5. **Impagliazzo-Wigderson (1997)**: If E requires 2^Omega(n) circuits, then BPP=P
6. **OWF implies PRG**: Blum-Micali-Yao construction via hardcore bits + stretch composition

### Core Algorithms

- **LCG**: X_{n+1} = (aX_n + c) mod m
- **MT19937**: Mersenne Twister, period 2^19937-1
- **Xorshift128**: Marsaglia, period 2^128-1
- **ChaCha20**: Bernstein, 20-round cryptographic stream cipher
- **NIST SP 800-22**: Frequency, Block, Runs, Longest Run, Serial, Entropy, Cumulative Sums
- **Berlekamp-Massey**: Linear complexity via LFSR synthesis over GF(2)
- **Floyd Algorithm**: Cycle detection in PRG state space

### Canonical Problems

- PRG distinguisher problem
- Next-bit prediction problem
- Hardcore predicate extraction
- OWF inversion / discrete log / factoring
- BPP derandomization via PRG

### Course Alignment

| School | Course | Topics |
|--------|--------|--------|
| MIT | 6.841 Adv Complexity | Yao Theorem, GL hardcore, IW97 |
| Stanford | CS254 Computational Complexity | PRG definitions, hybrid argument, hardcore predicates |
| Berkeley | CS278 Computational Complexity | BPP derandomization, NW generator |
| CMU | 15-855 Graduate Complexity | PRG constructions, LCG/MT analysis |
| Princeton | COS 522 Computational Complexity | Cryptographic PRGs, next-bit tests |
| Caltech | CS 151 Complexity Theory | Statistical tests, distinguisher bounds |
| Cambridge | Part III Advanced Complexity | OWF->PRG, hardcore bit extraction |
| Oxford | Advanced Complexity Theory | NIST test suite, randomness testing |
| ETH | 263-4650 Advanced Complexity | Formal PRG model, NW framework |

### File Structure

```
mini-pseudorandom-generators/
|-- Makefile                 # make test passes
|-- README.md                # This file
|-- include/                 # 6 header files (191 lines)
|   |-- prg_core.h           # Core PRG types and API
|   |-- prg_constructions.h  # LCG, MT19937, ChaCha20, DRBGs
|   |-- prg_constructs.h     # Lightweight state machine variants
|   |-- hardcore_predicates.h # GL, BM, BBS predicates + OWF
|   |-- next_bit_tests.h     # NIST test suite (10 tests)
|   |-- nist_tests.h         # Extended NIST tests + statistical tools
|-- src/                     # 6 implementation files (3411 lines)
|   |-- prg_core.c           # Core theorems: Yao, Hybrid, IW97, NW (569 lines)
|   |-- prg_constructions.c  # PRG implementations: LCG, MT, ChaCha20, DRBG (600 lines)
|   |-- hardcore_predicates.c # Hardcore predicates + OWF + analysis (662 lines)
|   |-- next_bit_tests.c     # NIST test implementations (646 lines)
|   |-- nist_tests.c         # Extended NIST tests (515 lines)
|   |-- prg_constructs.c     # State machine variants (349 lines)
|-- tests/
|   |-- test.c               # 54 comprehensive tests (683 lines)
|-- examples/
|   |-- example1.c           # Yao theorem + PRG security analysis
|   |-- example2.c           # PRG constructions + statistical testing
|   |-- example3.c           # Cryptographic applications + NIST suite
|-- docs/
    |-- knowledge-graph.md
    |-- coverage-report.md
    |-- gap-report.md
    |-- course-alignment.md
    |-- course-tree.md
```