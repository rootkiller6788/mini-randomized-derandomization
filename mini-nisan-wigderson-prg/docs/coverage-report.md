# Coverage Report — mini-nisan-wigderson-prg

## Knowledge Layer Assessment

| Level | Name | Status | Score | Evidence |
|-------|------|--------|-------|----------|
| **L1** | Definitions | **COMPLETE** | 2 | 10 core struct/typedef definitions across all headers |
| **L2** | Core Concepts | **COMPLETE** | 2 | 8 concepts: hard/random tradeoff, derandomization, indistinguishability, stretch, amplification |
| **L3** | Math Structures | **COMPLETE** | 2 | 10 structures: set systems, GF(q), DAG circuits, truth tables, Fourier analysis, expanders |
| **L4** | Fundamental Laws | **COMPLETE** | 2 | 10 theorems: NW Theorem, Yao XOR, Shannon, Lupanov, Impagliazzo, Hastad, RR Barrier, GL Theorem |
| **L5** | Algorithms/Methods | **COMPLETE** | 2 | 15 algorithms: RS design, random design, PRG eval, hybrid arg, distinguisher, BPP sim, statistical tests, KDF |
| **L6** | Canonical Problems | **COMPLETE** | 2 | 6 problems: PARITY, MAJORITY, MOD-p, Inner Product, PRG Distinguishing, Circuit LB |
| **L7** | Applications | **COMPLETE** | 2 | 5 apps: BPP derandomization, KDF, randomized sim, PRG optimization, stat. testing |
| **L8** | Advanced Topics | **COMPLETE** | 2 | 7 topics: Switching Lemma, Natural Proofs, Expanders, RR threshold, AC0 LB, Spectral, W2A |
| **L9** | Research Frontiers | **Partial** | 1 | 5 topics documented in knowledge graph (no implementation required per SKILL.md) |

**Total Score: 17/18**

## Code Metrics

| Metric | Value | Requirement |
|--------|-------|-------------|
| include/ + src/ lines | ~4,935 | ≥ 3,000 ✓ |
| include/*.h files | 5 | ≥ 4 ✓ |
| src/*.c files | 5 | ≥ 4 ✓ |
| src/*.lean files | 1 | ≥ 1 ✓ |
| tests/test.c | 814 lines | exists ✓ |
| examples/*.c files | 4 | ≥ 3 ✓ |

## Lean Formalization

- File: `src/nw_theorem.lean` (230 lines)
- Contains: NWDesign, HardFunction, NWPRG structures
- Theorem statements: Yao XOR Lemma, NW Theorem, Natural Proofs Barrier
- Uses `Bit` inductive type, `BitVec`, `TruthTable` definitions

## Test Coverage

- 37 test cases across L1-L8
- Covers: design creation/verification, PRG eval, XOR lemma, circuit construction, statistical tests, distinguisher advantage
- All tests pass with `make test`