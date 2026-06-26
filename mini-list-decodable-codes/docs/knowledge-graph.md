# Knowledge Graph — mini-list-decodable-codes

## L1: Core Definitions (Complete)
- Code parameters [n,k,d]_q
- Codeword, codeword list, Hamming distance
- List-decodable code: definition and list-size bound
- Reed-Solomon code: evaluation-point construction
- Linear code: generator matrix, parity-check matrix
- Cyclic code: generator polynomial
- LDPC code: sparse parity-check matrix
- Finite field GF(q): prime and extension fields
- Folded RS code: folding parameter, block alphabet
- Reed-Muller code RM(r,m): parameters n=2^m, k=ΣC(m,i), d=2^{m-r}
- Hadamard code H_k: n=2^k, k, d=2^{k-1}
- Univariate and bivariate polynomials over GF(q)

## L2: Core Concepts (Complete)
- Encoding: message → codeword (polynomial eval / matrix mult)
- List-decoding radius (Johnson radius)
- Combinatorial vs algorithmic list-decoding
- Unique decoding limit (d/2) vs list-decoding beyond
- Code rate R = k/n, relative distance δ = d/n
- Polynomial list size (polynomially bounded in n)
- Johnson bound: e ≤ n(1-sqrt(1-d/n)) ⇒ L ≤ qn
- Capacity: δ < 1-R-ε implies polynomial list size
- Systematic encoding for linear codes
- Belief propagation for LDPC codes

## L3: Mathematical Structures (Complete)
- Finite field GF(p): modular arithmetic (add, sub, mul, inv, pow)
- Finite field GF(p^m): polynomial representation
- Primitive element and multiplicative group structure
- Polynomial ring GF(q)[x]: degree, addition, multiplication
- Euclidean algorithm for polynomials over GF(q)
- Lagrange interpolation over finite fields
- Chien search: root-finding in GF(q)
- Bivariate polynomials and (1,k-1)-weighted degree
- Hasse derivative: multiplicity constraints for GS algorithm
- Monomial ordering for Sudan interpolation
- Tanner graph for LDPC codes (variable/check nodes)

## L4: Fundamental Laws (Complete)
- Singleton bound: k ≤ n-d+1 (MDS codes achieve it)
- Hamming (sphere-packing) bound
- Plotkin bound: low-rate regime δ > (q-1)/q
- Gilbert-Varshamov bound: existence R ≥ 1-H_q(δ)
- Griesmer bound: Σ ⌈d/q^i⌉ ≤ n for linear codes
- Johnson bound: list-size bound for list decoding
- Zyablov-Pinsker capacity: δ_cap = 1-R
- MRRW (LP) bound: best known upper bound for binary codes
- Elias-Bassalygo bound
- Relaxed Singleton bound for list decoding
- List-recovery capacity bound
- Johnson bound for RS codes: e_J = n(1-sqrt{R})

## L5: Algorithms/Methods (Complete)
- Reed-Solomon encoding via Horner polynomial evaluation
- Welch-Berlekamp unique decoder (rational interpolation)
- Sudan list-decoding algorithm (1997):
  - Bivariate interpolation with (1,k-1)-weighted degree
  - Linear-algebraic root-finding
- Guruswami-Sudan multiplicity-based algorithm (1999):
  - Hasse derivative constraints
  - Improved decoding radius
- Goldreich-Levin algorithm (1989):
  - Local list-decoding of Hadamard codes
  - Pairwise independence + prefix guessing
  - Hard-core predicate extraction
- Gopalan-Klivans-Zuckerman list-decoding for RM codes
- Reed's majority-logic decoder for RM codes
- Belief propagation (sum-product) for LDPC codes
- Linear-algebraic folded RS decoder (Guruswami-Xing 2012)
- Lagrange interpolation over GF(q)
- Fast exponentiation (square-and-multiply) over GF(q)

## L6: Canonical Problems (Complete)
- Reed-Solomon code list-decoding (Sudan, GS)
- Folded RS code capacity-achieving list-decoding
- Hadamard code list-decoding (Goldreich-Levin)
- Reed-Muller code list-decoding (GKZ)
- Linear code bounded-distance and list decoding
- BCH code syndrome decoding
- LDPC belief-propagation decoding
- Unique decoding (Welch-Berlekamp) baseline
- Enumerating all codewords for small codes (ground truth)

## L7: Applications (Complete)
- Hardness amplification: (1-ε)-hard → (1/2+ε')-hard via codes
  (Sudan-Trevisan-Vadhan 2001, Trevisan 2003)
- Randomness extraction from list-decodable codes
  (Trevisan extractor paradigm)
- Soft-decision decoding (Koetter-Vardy 2003 framework)
- Hard-core predicate extraction (Goldreich-Levin)
  from one-way functions
- Derandomization via list-decodable codes

## L8: Advanced Topics (Complete)
- Capacity-achieving codes: folded RS proof structure
- Concatenated code construction (Forney 1966)
  outer list-decodable + inner small-alphabet code
- List recoverability: S_i input sets, output c_i∈S_i
- Deterministic Goldreich-Levin using ε-biased sets (Naor-Naor)
- Fourier analysis of Boolean functions (heavy coefficients)
- LP lower bounds for codes (MRRW)

## L9: Research Frontiers (Complete)
- Approximate (relaxed) list-decoding: `ld_relaxed_list_decode()`
  Guruswami 2003 — decoding beyond unique radius via relaxation
- Subspace evasive sets: `ses_construct_polynomial()`, `ses_verify_evasive()`
  Dvir-Lovett 2012 — construction and verification
- Quantum CSS codes: `qcss_create_from_classical()`, `qcss_list_decode_errors()`
  Calderbank-Shor-Steane construction from list-decodable classical codes
- Meta-complexity / MCSP: `mscp_encode_circuit()`, `mscp_list_decode_circuits()`
  Kabanets-Cai 2000 — circuit size minimization via list-decoding
- LTC + list-decoding pipeline: `ltc_test_linearity()`, `ltc_list_decode_pipeline()`
  Hadamard code as both locally testable (BLR) and list-decodable (GL)
