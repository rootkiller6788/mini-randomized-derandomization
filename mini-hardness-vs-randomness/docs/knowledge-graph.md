# Knowledge Graph ¡ª mini-hardness-vs-randomness

## L1: Definitions
- Hard function (¦Ä-hard against circuits of size S)
- Pseudorandom generator (PRG)
- Circuit lower bound certificate
- Hardness level (WEAK, MODERATE, STRONG, EXPONENTIAL)
- Nisan-Wigderson PRG
- Boolean circuit DAG
- Gate types (AND, OR, NOT, XOR, MAJORITY, MODp)
- BPP algorithm
- Hardness certificate
- Hardcore set
- Worst-case hardness
- Average-case hardness

## L2: Core Concepts
- Hardness vs Randomness paradigm (Nisan-Wigderson 1994)
- If hard functions exist ¡ú PRG exists ¡ú BPP = P
- Hardness amplification (worst-case ¡ú average-case)
- Derandomization via hardness
- Circuit complexity as hardness measure
- The Nisan-Wigderson construction
- Impagliazzo's five worlds
- Pseudorandomness from hardness
- BPP ? P/poly (Adleman)
- One-way functions from hardness

## L3: Mathematical Structures
- Boolean circuits as DAGs of gates
- Circuit families and uniformity
- Gate-level evaluation
- Circuit depth, size, fan-in/fan-out
- Truth tables (2^n entries for n-variable function)
- Hardness certificates (n, circuit_sz, hardness)
- PRG parameters (seed_len, output_len, stretch)
- NW combinatorial designs
- Hardcore sets (density-based)
- WTA conversion parameters

## L4: Fundamental Laws/Theorems
- Shannon's counting lower bound (1949)
- Lupanov upper bound (1958)
- H?stad's switching lemma (1986) - PARITY ? AC0
- Smolensky's algebraic lower bound (1987) - MODp ? AC0[q]
- Razborov's monotone circuit lower bound (1985)
- Alon-Boppana refinement (1987)
- Furst-Saxe-Sipser / Yao (1984/1985) - MAJORITY ? AC0
- Nisan-Wigderson Theorem (1994): Hardness ¡ú PRG ¡ú Derandomization
- Impagliazzo-Wigderson Theorem (1997): BPP = P under exponential hardness
- Impagliazzo's Hardcore Lemma (1995)
- Yao's XOR Lemma (1982)
- Goldreich-Levin Theorem (1989)
- Adleman's Theorem: BPP ? P/poly (1978)

## L5: Algorithms/Methods
- Nisan-Wigderson PRG construction
- Yao's XOR lemma amplification
- Impagliazzo hardcore set construction
- Conditional expectations derandomization
- Adleman derandomization (advice string)
- NW combinatorial design construction
- PRG stretching computation
- Hardness certification from circuit lower bounds
- Error probability analysis
- Deterministic simulation of BPP

## L6: Canonical Problems
- PARITY function (AC0 lower bound)
- MAJORITY function (AC0 lower bound)
- CLIQUE problem (monotone circuit lower bound)
- MODp computation (algebraic lower bound)
- PRG construction from hard predicates
- BPP derandomization
- Hardness amplification pipeline
- Circuit evaluation
- Shannon bound computation
- Formula complexity estimation

## L7: Applications
- One-way functions from hardness (cryptography)
- Deterministic primality testing (AKS vs Miller-Rabin)
- Cryptographic key size selection
- GPS pseudorandom codes (Gold codes)
- Derandomization in practice
- NASA Monte Carlo simulation derandomization
- Climate model ensemble methods
- Nuclear reactor safety simulation (IAEA standards)
- Deterministic scenario generation (autonomous driving)
- ISO 26262 functional safety verification

## L8: Advanced Topics
- Impagliazzo's five worlds framework
- Exponential Complexity Assumption (ECA)
- Non-black-box derandomization
- Pseudo-deterministic algorithms
- ETH ¡ú IW theorem gap analysis
- NW design combinatorial construction
- List decoding hardness
- Hardcore predicate extraction

## L9: Research Frontiers
- Meta-complexity connections
- Natural proofs barrier relevance
- Quantum derandomization
- Fine-grained hardness-vs-randomness
- GCT and circuit lower bounds