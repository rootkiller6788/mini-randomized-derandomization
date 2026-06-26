# Knowledge Graph — mini-expander-graphs-construction

## L1: Definitions
- ExpanderGraph (d-regular graph, adjacency list, spectral data)
- GraphVector (real-valued function on vertices)
- RamanujanParams (p, q primes for LPS construction)
- ZigZagParams / ZigZagResult (product graph data)
- RandomWalk (graph + step count)
- Spectral expander definition: λ₂ ≤ (1-γ)·d
- Vertex expansion, edge expansion (Cheeger constant)
- Ramanujan graph: λ₂ ≤ 2√(d-1)

## L2: Core Concepts
- Spectral gap γ = d - λ₂
- Mixing time τ(ε) of lazy random walk
- Conductance / Cheeger constant h(G)
- Vertex expansion vs edge expansion
- Explicit vs probabilistic constructions
- Lazy random walk (stay with prob 1/2)
- Total variation distance from uniform
- Collision probability, entropy rate

## L3: Mathematical Structures
- Adjacency matrix A
- Laplacian matrix L = D - A
- Normalized Laplacian L_norm = I - A/d
- Rotation map Rot: [n]×[d] → [n]×[d]
- Cayley graphs over Z_n
- Tensor product, graph power
- Zig-zag product: G(z)H on n_G·d_G vertices
- Replacement product: G(r)H
- PGL(2,q) group for LPS construction
- Quadratic residues, Legendre symbol

## L4: Fundamental Theorems
- Alon-Boppana bound: liminf λ₂ ≥ 2√(d-1) for any d-regular family
- Nilli (finite-n) bound: quantitative version of Alon-Boppana
- Cheeger inequality: γ/2 ≤ h(G) ≤ √(2dγ) for regular graphs
- Expander Mixing Lemma (Alon-Chung 1988): |E(S,T) - (d/n)|S||T|| ≤ λ√(|S||T|)
- Tanner's inequality / Hoffman bound: α(G) ≤ n·λ₂/(d+λ₂)
- RVW zig-zag spectral bound: λ(G(z)H) ≤ λ_G + λ_H + λ_H²
- Expander Chernoff bound (Gillman 1998)
- Completeness of LPS construction: λ₂ ≤ 2√(d-1) for p,q ≡ 1 mod 4

## L5: Algorithms/Methods
- Power iteration (deflated) for λ₂ computation
- Jacobi eigenvalue algorithm (full spectrum)
- Margulis-Gabber-Galil explicit construction (O(n))
- LPS Ramanujan construction via quaternion algebra
- Zig-zag product construction (O(n·d·d_H²))
- Replacement product construction
- Random d-regular via configuration model
- Cayley graph construction on Z_n
- Lazy random walk simulation
- Spectral partitioning via Fiedler vector
- BFS connectivity and distance computation

## L6: Canonical Problems
- USTCON: undirected s-t connectivity on expanders
- Error amplification for RP/BPP algorithms
- Derandomized sampling on expanders
- Expander sorting network (AKS framework)
- Expander code encoding/decoding (Sipser-Spielman)
- PRG generation via expander walks

## L7: Applications
- Error reduction for randomized algorithms (saves random bits)
- Derandomized sampling with statistical guarantees
- Expander error-correcting codes (linear-time encoding/decoding)
- Pseudorandom generators from expander walks
- Extractor construction based on expanders
- Load balancing on expander interconnection networks
- Deterministic log-space USTCON (SL = L, Reingold 2008)

## L8: Advanced Topics
- Heat kernel and continuous-time diffusion on expanders
- Lossless expanders (neighbor set nearly size d·|S|)
- Full eigenvalue spectrum via Jacobi method
- Unique neighbor expanders (conceptual)
- High-dimensional expanders (conceptual, documented)

## L9: Research Frontiers
- Coboundary expanders and agreement testing
- High-dimensional expanders for PCP improvement
- Quantum expanders and applications
- Explicit constructions of lossless expanders