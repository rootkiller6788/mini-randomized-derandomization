# Coverage Report — mini-randomness-extractors

## Summary

| Level | Status | Score | Evidence |
|-------|--------|-------|----------|
| L1 Definitions | **COMPLETE** | 2 | 10+ typedef structs, Lean definitions |
| L2 Core Concepts | **COMPLETE** | 2 | 10 entropy/distance measures implemented |
| L3 Math Structures | **COMPLETE** | 2 | GF(2), expander graphs, hash families, designs |
| L4 Fundamental Laws | **COMPLETE** | 2 | LHL, RTD, NZ, Chernoff, Pinsker all implemented |
| L5 Algorithms | **COMPLETE** | 2 | 12 extractor constructions implemented |
| L6 Canonical Problems | **COMPLETE** | 2 | 8 problems with examples |
| L7 Applications | **COMPLETE** | 2 | 6 applications: derandomization, crypto, QKD |
| L8 Advanced Topics | **PARTIAL** | 1 | 6 topics: Trevisan, Bourgain, quantum, network |
| L9 Research Frontiers | **PARTIAL** | 1 | Documented: quantum, network, non-malleable |

**Total: 16/18 — COMPLETE**

## Level Details

### L1: Complete
- 10+ typedef struct definitions across 4 headers
- All core types instantiated in C and formalized in Lean
- `Extractor`, `WeakSource`, `SeededExt`, `TwoSourceExt`, `Disperser`, `ObliviousDisp`, `UniversalHash`, `ExtractorGraph`

### L2: Complete
- Min-entropy, Shannon, Renyi (variable alpha) entropy
- Statistical distance, Hellinger, KL divergence
- Entropy deficit
- All have tests with mathematical assertions

### L3: Complete
- GF(2) field: add (XOR), multiply (polynomial mod irreducible)
- Expander graphs: spectral expansion, neighbor fraction
- Carter-Wegman universal hash (prime field)
- Pairwise independent hash (GF(2) matrix-vector)
- Combinatorial designs (polynomial method)

### L4: Complete
- Leftover Hash Lemma: error bound computation + Lean statement
- RTD lower bound: output length validation
- NZ seed bound: seed length computation
- Chernoff: min-entropy verification with confidence interval
- Pinsker, Gibbs, triangle inequality tested

### L5: Complete
- 5 seeded extractors: LHL, Trevisan, Raz, Ta-Shma, block-source
- 2 two-source extractors: Chor-Goldreich, Bourgain
- Extractor graph construction + verification
- NW PRG from hard predicates
- Carter-Wegman and pairwise independent hashing
- Disperser construction (polynomial + expander)

### L6: Complete
- Extract from biased coin, estimate distance, verify entropy
- Build custom extractors, 2-source extraction
- Derandomize BPP, key derivation, privacy amplification
- All demonstrated in 3 examples + test suite

### L7: Complete (3+ applications)
1. BPP derandomization (NW/IW paradigm)
2. Cryptographic key derivation from weak secrets
3. Privacy amplification (BBR88 protocol)
4. Explicit extractor simulation with statistical validation
5. Network extractor protocol
6. Statistical test battery

### L8: Partial
1. Trevisan extractor (optimal seed, NW+designs) ✓
2. Bourgain 2-source extractor ✓
3. Extractor graphs with spectral bounds ✓
4. Multi-source fault tolerance ✓
5. Quantum-proof extractor bounds ✓
6. Network extractor protocols ✓

### L9: Partial (documented)
- Quantum-proof extraction bounds
- Network extractor research (KRV 2009)
- Non-malleable extractors
- Multi-source extraction
