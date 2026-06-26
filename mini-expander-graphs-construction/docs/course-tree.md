# Prerequisite Tree — mini-expander-graphs-construction

## Dependencies

### Immediate Prerequisites
- **mini-complexity-foundations**: Basic complexity classes (P, NP, RP, BPP, SL, L)
- **mini-boolean-circuits-model**: Circuit complexity background for expander-based derandomization connections
- **Linear algebra**: Eigenvalues, eigenvectors, matrix operations
- **Graph theory**: Basic definitions, adjacency matrices, degrees
- **Number theory**: Modular arithmetic, quadratic residues (for LPS construction)
- **Probability theory**: Random walks, Markov chains, mixing times

### Knowledge Flow
```
Linear Algebra ──► Spectral Graph Theory ──► Expander Analysis
       │                                            │
Graph Theory  ──► Regular Graphs ──► Explicit Constructions
       │                                            │
Probability  ──► Random Walks ──► Mixing Time, Connectivity
                                                    │
Complexity Theory ──► Derandomization ◄─────────────┘
                     SL = L (Reingold 2008)
                     Error Reduction for RP/BPP
```

### Downward Dependencies
This module is a prerequisite for:
- **mini-pcp-theorem**: Expanders are a key ingredient in PCP constructions
- **mini-ip-pspace**: Interactive proofs use expander-based sampling
- **mini-communication-complexity**: Expander-based protocols
- **mini-approximation-algorithms**: Hardness of approximation via expander graphs

### Research Frontier Connections (L9)
- Coboundary expanders → PCP Theorem improvements
- High-dimensional expanders → Agreement testing, property testing
- Quantum expanders → Quantum error correction, quantum derandomization
- Explicit lossless expanders → Coding theory, extractors