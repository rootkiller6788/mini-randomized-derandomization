# mini-randomness-extractors

## Module Status: COMPLETE

| Level | Status |
|-------|--------|
| L1 Definitions | Complete |
| L2 Core Concepts | Complete |
| L3 Mathematical Structures | Complete |
| L4 Fundamental Laws | Complete |
| L5 Algorithms/Methods | Complete |
| L6 Canonical Problems | Complete |
| L7 Applications | Complete (6 applications) |
| L8 Advanced Topics | Partial (6 topics) |
| L9 Research Frontiers | Partial (documented) |

**Total: 16/18 (COMPLETE)** · **Line count (include/ + src/): 3330+ (qualified)**

---

## Core Definitions

A **randomness extractor** is a function `Ext : {0,1}^n × {0,1}^d → {0,1}^m` such that for any random variable X over {0,1}^n with min-entropy H∞(X) ≥ k, the output Ext(X, U_d) is ε-close to the uniform distribution U_m in statistical distance.

**Types:**
- `Extractor` — (k,ε)-extractor parameters
- `WeakSource` — k-source with min-entropy guarantee
- `SeededExt` — extractor with embedded function pointer
- `TwoSourceExt` — independent-source extractor
- `Disperser` — weaker relaxation: hits almost entire range
- `ObliviousDisp` — zero-seed disperser
- `UniversalHash` — Carter-Wegman hash family
- `ExtractorGraph` — bipartite expander representation

---

## Core Theorems

### Leftover Hash Lemma (Impagliazzo, Levin, Luby 1989)
If H is a universal hash family from {0,1}^n to {0,1}^m and H∞(X) ≥ k, then:
```
Δ((H, H(X)), (H, U_m)) ≤ 2^{-(k-m)/2}
```

### RTD Bound (Radhakrishnan, Ta-Shma 2000)
Any (k,ε)-extractor must satisfy:
```
m ≤ k + 2·log₂(ε) + 2
```

### NZ Seed Bound (Nisan, Zuckerman 1996)
```
d ≥ max(1, log₂(n-k) + 2·log₂(1/ε) - O(1))
```

### Chernoff Inequality
For N i.i.d. Bernoulli trials with mean μ:
```
Pr[|X̄ - μ| ≥ ε] ≤ 2·exp(-2Nε²)
```

### Chor-Goldreich (1988)
If X,Y are independent n-bit sources each with min-entropy > n/2:
```
Δ(⟨X,Y⟩ mod 2, U₁) ≤ 2^{-(H∞(X)+H∞(Y)-n-1)/2}
```

---

## Core Algorithms

| Algorithm | Construction | Seed Length | Output | Reference |
|-----------|-------------|-------------|--------|-----------|
| LHL | Universal hashing | n | k - 2log(1/ε) | ILL 1989 |
| Trevisan | NW PRG + designs | O(log n) | k^Ω(1) | Trevisan 2001 |
| Raz | Block-and-condense | O(log n) | ~k | Raz 2005 |
| Ta-Shma | Code-based (RS) | O(log n) | k - O(1) | Ta-Shma 1996 |
| Chor-Goldreich | Inner product GF(2) | 0 (2-source) | 1 bit (basic) | CG 1988 |
| Bourgain | Field multiplication | 0 (2-source) | Ω(n) bits | Bourgain 2005 |

---

## Canonical Problems

1. Extract unbiased bits from a biased coin source
2. Estimate statistical distance of a distribution to uniform
3. Verify min-entropy bound from empirical samples
4. Build a (k,ε)-extractor for given source parameters
5. Perform 2-source extraction from independent weak sources
6. Derandomize a BPP algorithm using seed enumeration
7. Derive a cryptographic key from a Diffie-Hellman shared secret
8. Execute privacy amplification after quantum key distribution

---

## Nine-School Course Mapping

| School | Course | Topic Alignment |
|--------|--------|----------------|
| MIT | 6.841 Advanced Complexity | Extractors, dispersers, pseudorandomness |
| Stanford | CS254 Computational Complexity | Derandomization, NW PRG |
| Berkeley | CS278 Computational Complexity | PCP, extractors, samplers |
| CMU | 15-855 Graduate Complexity | Full extractor theory |
| Princeton | COS 522/551 | Cryptographic applications |
| Caltech | CS 151/154 | Extractor lower bounds, NW |
| Cambridge | Part II/III Complexity | Core definitions, 2-source |
| Oxford | Advanced Complexity Theory | Expander-based extractors |
| ETH | 263-4650 Advanced Complexity | Pseudorandomness, formal proofs |

---

## Project Structure

```
mini-randomness-extractors/
├── Makefile                    # make test (compiles + runs)
├── README.md                   # This file
├── PROMPT.md                   # Task specification
├── include/                    # 4 header files (859 lines)
│   ├── extractor_core.h        # Core types + entropy measures
│   ├── extractor_constructions.h # Extractor constructions + hashing
│   ├── extractor_applications.h  # Applications (BPP, crypto, QKD)
│   └── disperser_core.h        # Disperser definitions
├── src/                        # 5 source files (2471+ lines)
│   ├── extractor_core.c        # Entropy, distance, GF(2), extractor eval
│   ├── extractor_constructions.c # LHL, Trevisan, Raz, Ta-Shma, 2-source
│   ├── extractor_applications.c  # BPP derand, key-deriv, privacy amp
│   ├── disperser_core.c        # Disperser constructions + verification
│   └── extractor_theorems.lean # Lean 4 formalization (148 lines)
├── tests/
│   └── test.c                  # Comprehensive test suite (569 lines, 0 FAIL)
├── examples/
│   ├── example1.c              # Entropy & statistical distance demo
│   ├── example2.c              # Extractor constructions demo
│   └── example3.c              # Applications demo
├── docs/
│   ├── knowledge-graph.md      # L1-L9 knowledge coverage
│   ├── coverage-report.md      # Status per level
│   ├── gap-report.md           # Missing items + priority
│   ├── course-alignment.md     # Nine-school course mapping
│   └── course-tree.md          # Prerequisite dependencies
├── demos/                      # Visualization/演示
└── benches/                    # Performance benchmarks
```

## Building and Testing

```bash
make          # Build and run all tests (0 failures expected)
make clean    # Remove build artifacts
make check-lines  # Show line count for include/ + src/
```

## Key References

1. Impagliazzo, Levin, Luby (1989). "Pseudo-random Generation from One-way Functions." STOC 1989.
2. Nisan, Zuckerman (1996). "Randomness is Linear in Space." JCSS.
3. Trevisan (2001). "Extractors and Pseudorandom Generators." JACM.
4. Radhakrishnan, Ta-Shma (2000). "Bounds for Dispersers, Extractors, and Depth-Two Superconcentrators." SIAM J. Disc. Math.
5. Chor, Goldreich (1988). "Unbiased Bits from Sources of Weak Randomness." SIAM J. Comput.
6. Bourgain (2005). "More on the Sum-Product Phenomenon..." GAFA.
7. Bennett, Brassard, Robert (1988). "Privacy Amplification by Public Discussion." SIAM J. Comput.
8. Vadhan (2012). "Pseudorandomness." Foundations and Trends in TCS.
9. Arora, Barak (2009). "Computational Complexity: A Modern Approach." Chapters 20-21.
10. Goldreich (2008). "Computational Complexity: A Conceptual Perspective."
