# mini-hardness-vs-randomness

## Module Status: COMPLETE ?

- **L1-L7: Complete**
- **L8: Partial** (8 advanced topics, mostly implemented)
- **L9: Partial** (research frontiers documented)

**include/ + src/ = 3234 lines** | **97 tests, 0 failures** | **7 C sources + 1 Lean formalization**

## Core Paradigm

If computationally hard functions exist (provable circuit lower bounds), then pseudorandom generators exist, which can derandomize all of BPP ¡ª proving P = BPP. This is the Nisan-Wigderson / Impagliazzo-Wigderson paradigm.

### Core Theorem Chain
1. **Shannon (1949)**: Almost all functions require exponential-size circuits
2. **Yao (1982)**: XOR Lemma ¡ª worst-case hardness amplifies to average-case
3. **H?stad (1986)**: PARITY requires exponential AC0 circuits
4. **Razborov (1985)**: CLIQUE requires large monotone circuits
5. **Smolensky (1987)**: MODp requires exponential AC0[q] circuits
6. **Nisan-Wigderson (1994)**: Hardness implies pseudorandom generators
7. **Impagliazzo-Wigderson (1997)**: P = BPP if E requires exponential circuits

## Key Definitions (L1)

| Term | Definition |
|------|------------|
| Hard function | f: {0,1}^n ¡ú {0,1} requiring circuits of size S to compute |
| ¦Å-hardness | Every circuit of size S fails on ¡Ý ¦Å fraction of inputs |
| PRG | G: {0,1}^k ¡ú {0,1}^m, m > k, output computationally indistinguishable from random |
| Circuit lower bound | Proof that a function requires ¡Ý S gates in a given circuit class |
| Hardness certificate | (n, circuit_sz, hardness) triple certifying a lower bound |
| Hardcore set | Subset of inputs where function is average-case hard |
| Five worlds | Algorithmica, Heuristica, Pessiland, Minicrypt, Cryptomania |

## Core Theorems (L4)

| Theorem | Statement | Implementation |
|---------|-----------|----------------|
| Shannon (1949) | Most n-var Boolean functions require ¡Ý 2^n/(3n) gates | `clb_shannon_counting()` |
| H?stad (1986) | PARITY requires ¡Ý 2^{¦¸(n^{1/(d-1)})} size AC0 circuits | `clb_hastad_switching()` |
| Razborov (1985) | CLIQUE requires ¡Ý 2^{¦¸(k^{1/3})} monotone gates | `clb_razborov_clique()` |
| Smolensky (1987) | MODp requires exponential AC0[q] size (p¡Ùq primes) | `clb_smolensky_modp()` |
| NW Theorem (1994) | Hardness ¡ú PRG ¡ú Derandomization | `hr_implies_prg()` |
| IW Theorem (1997) | E requires exponential circuits ¡ú P = BPP | `iw_theorem_conditions_met()` |
| Yao XOR Lemma (1982) | k-fold XOR amplifies hardness exponentially | `wta_yao_xor_lemma()` |
| Impagliazzo Hardcore (1995) | ¦Ä-hard function has hardcore set of density ¦Ä | `wta_impagliazzo_hardcore()` |

## Core Algorithms (L5)

1. **NW PRG Construction**: Hard function + combinatorial design ¡ú pseudorandom generator
2. **Yao XOR Amplification**: Iterate XOR to convert worst-case to average-case hardness
3. **Hardcore Set Construction**: Identify dense subset where function is unpredictable
4. **Conditional Expectations**: Deterministically find good random string
5. **Adleman Derandomization**: Fixed advice string works for all inputs (BPP ? P/poly)

## Canonical Problems (L6)

- **PARITY**: Computing XOR of n bits ¡ª exponential AC0 lower bound
- **MAJORITY**: Threshold function ¡ª exponential AC0 lower bound
- **CLIQUE**: Graph clique detection ¡ª exponential monotone lower bound
- **MODp**: Modular counting ¡ª exponential AC0[q] lower bound for p¡Ùq

## Applications (L7)

| Application | Domain | Key Reference |
|-------------|--------|---------------|
| One-way functions from hardness | Cryptography | HILL (1999) |
| Deterministic primality testing | Number Theory | AKS (2002) |
| Cryptographic key sizing | Security | NIST SP 800-57 |
| GPS Gold codes | Signal Processing | IS-GPS-200 |
| NASA Monte Carlo derandomization | Simulation | NASA CR-188243 |
| Climate ensemble methods | Climate Science | IPCC AR6 |
| Nuclear safety analysis | Nuclear Engineering | IAEA 2015 |
| ISO 26262 verification | Automotive Safety | ISO 26262:2018 |

## Nine-School Course Mapping

| School | Course | Topics Mapped |
|--------|--------|---------------|
| MIT | 6.841 Adv Complexity | Hardness vs Randomness (Weeks 7-9) |
| Stanford | CS254 Comp Complexity | Circuit lower bounds, derandomization |
| Berkeley | CS278 Comp Complexity | PCP-adjacent derandomization |
| CMU | 15-855 Grad Complexity | Circuit complexity, PRG constructions |
| Princeton | COS 522 Comp Complexity | Cryptographic applications |
| Caltech | CS 151 Complexity Theory | Hardness vs Randomness (Week 8) |
| Cambridge | Part III Adv Complexity | Circuit complexity + derandomization |
| Oxford | Adv Complexity Theory | Randomness and Derandomization |
| ETH | 263-4650 Adv Complexity | Logical approaches, five worlds |

## File Structure

```
mini-hardness-vs-randomness/
©À©¤©¤ Makefile                         # make test ¡ª 97 tests, all pass
©À©¤©¤ README.md                        # This file
©À©¤©¤ include/
©¦   ©À©¤©¤ circuit_lower_bounds.h       # Boolean circuits, lower bound API
©¦   ©À©¤©¤ hardness_randomness.h        # Hardness ¡ú PRG connection
©¦   ©À©¤©¤ derandomization_via_hardness.h # BPP derandomization API
©¦   ©¸©¤©¤ worst_to_average.h           # Hardness amplification API
©À©¤©¤ src/
©¦   ©À©¤©¤ circuit_lower_bounds.c       # Circuit eval, Shannon, H?stad, Razborov
©¦   ©À©¤©¤ hardness_randomness.c        # HR connection, NW design
©¦   ©À©¤©¤ derandomization_via_hardness.c # BPP derandomization, Adleman
©¦   ©À©¤©¤ worst_to_average.c           # Yao XOR, Impagliazzo hardcore
©¦   ©À©¤©¤ iw_theorem.c                 # IW theorem, five worlds
©¦   ©À©¤©¤ applications.c               # Crypto, GPS, NASA, nuclear, ISO
©¦   ©À©¤©¤ circuit_simulation.c         # Circuit construction, analysis
©¦   ©¸©¤©¤ core.lean                    # Lean 4 formalization (15+ theorems)
©À©¤©¤ tests/
©¦   ©¸©¤©¤ test.c                       # 97 assert-based tests
©À©¤©¤ examples/
©¦   ©À©¤©¤ example_circuit_lb.c         # Circuit lower bound demo
©¦   ©À©¤©¤ example_derandomization.c    # Derandomization pipeline demo
©¦   ©¸©¤©¤ example_hardness_amplify.c   # Hardness amplification demo
©À©¤©¤ benches/
©¦   ©¸©¤©¤ bench_circuit_lb.c           # Performance benchmarks
©À©¤©¤ demos/
©¦   ©¸©¤©¤ demo.c                       # Interactive demonstration
©¸©¤©¤ docs/
    ©À©¤©¤ knowledge-graph.md           # L1-L9 knowledge coverage table
    ©À©¤©¤ coverage-report.md           # Detailed coverage analysis
    ©À©¤©¤ gap-report.md                # Remaining gaps
    ©À©¤©¤ course-alignment.md          # Nine-school course mapping
    ©¸©¤©¤ course-tree.md               # Prerequisite dependency tree
```

## Key References

- Arora & Barak (2009): *Computational Complexity: A Modern Approach*, Chapters 6, 14, 19, 20
- Nisan & Wigderson (1994): *Hardness vs Randomness*, JCSS 49(2)
- Impagliazzo & Wigderson (1997): *P = BPP if E requires exponential circuits*, JCSS
- Impagliazzo (1995): *A personal view of average-case complexity*, CCC
- Yao (1982): *Theory and applications of trapdoor functions*, FOCS
- H?stad (1986): *Almost optimal lower bounds for small depth circuits*, STOC
- Razborov (1985): *Lower bounds for the monotone complexity of Boolean functions*
- Goldreich & Levin (1989): *Hard-core predicate for any one-way function*, STOC
- Vollmer (1999): *Introduction to Circuit Complexity*
- Jukna (2012): *Boolean Function Complexity*