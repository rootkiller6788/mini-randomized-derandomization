# Course Tree — mini-randomness-extractors

## Prerequisites

```
Probability Theory (basic)
  └── Entropy measures (H_inf, H, D_KL)
      └── min-entropy, statistical distance

Complexity Theory (basic)
  └── P, BPP, RP definitions
      └── BPP derandomization

Linear Algebra (basic)
  └── GF(2) field operations

Graph Theory (basic)
  └── Expander graphs, spectral expansion
```

## Dependency Graph within Module

```
L1: Definitions
 ├── L2: Entropy measures
 │    ├── L4: LHL, RTD, NZ bounds
 │    │    └── L5: LHL extractor
 │    └── L3: GF(2), expanders
 │         ├── L5: Trevisan, Raz, Ta-Shma
 │         └── L6: 2-source extractors
 ├── L2: Statistical distance
 │    └── L4: Pinsker, triangle inequality
 └── L1: Disperser
      ├── L5: Polynomial disperser
      └── L6: Expander-based disperser

L7: Applications
 ├── L5: extractors + L4: bounds → key derivation
 ├── L5: NW PRG → BPP derandomization
 └── L4: LHL → privacy amplification

L8: Advanced Topics
 ├── L5: Trevisan → optimal seed extractors
 ├── L5: 2-source → Bourgain sum-product
 └── L3: expanders → extractor graphs

L9: Research Frontiers
 ├── L8: quantum → quantum-proof extraction
 ├── L8: network → multi-party protocols
 └── (External) non-malleability, seedless extraction
```

## Course Mapping (External Dependencies)

```
CS Theory Core
 ├── Automata Theory (Sipser)
 │    └── Turing machines → L1 Extractor definition
 ├── Computational Complexity (Arora-Barak)
 │    ├── Chapter 20: Derandomization → L7 BPP
 │    └── Chapter 21: Extractors → L4+L5
 ├── Probability & Computing (Mitzenmacher-Upfal)
 │    ├── Chernoff bounds → L4 verification
 │    └── Entropy → L2 measures
 └── Cryptography (Katz-Lindell)
      └── Hash functions → L5 universal hashing
```

## Next Modules

This module is a dependency for:
- `mini-approximation-algorithms/` (hardness of approximation)
- `mini-pcp-theorem/` (PCP uses randomness + verification)
- `mini-ip-pspace/` (interactive proofs use randomness)
- `mini-communication-complexity/` (randomized protocols)
