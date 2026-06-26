# Knowledge Graph — BPP/RP/ZPP Classes

## L1: Definitions (Complete ✅)
- RandomSource: randomness abstraction with seed, state, counter
- ProbabilisticTM: probabilistic Turing machine with random tape
- PTMOutcome: correct accept/reject, false accept/reject
- RandomizedClass: BPP, RP, co-RP, ZPP, P, UNKNOWN
- RandomizedDecision: confidence-bounded decision result
- GateType: enumerated Boolean gate types
- CircuitGate: single gate with type, inputs, evaluation
- BooleanCircuit: DAG circuit model with random bit inputs
- ProbabilisticCircuitFamily: indexed family {C_n}
- RandomizedReduction: BPP/RP reduction between languages
- RandomizedLanguage: language with decider and class

## L2: Core Concepts (Complete ✅)
- BPP amplification via majority vote (Chernoff bound)
- RP amplification via OR of trials
- ZPP simulation with early termination
- Error probability bounds (1/3 base, 2^{-k} amplified)
- PTM simulation and outcome classification
- Monte Carlo acceptance probability estimation
- Randomized reduction composition and verification

## L3: Mathematical Structures (Complete ✅)
- Probabilistic TM formalization (Q,Σ,Γ,δ with random tape)
- Boolean circuit DAG model with random inputs
- Probability distributions on {0,1}^m
- Pairwise independent hash families (Carter-Wegman)
- k-wise independent hash families (polynomials over GF(p))
- Chernoff bound: multiplicative form, KL-divergence
- Hoeffding bound: additive form
- Azuma-Hoeffding: martingale concentration
- McDiarmid's inequality: bounded differences

## L4: Fundamental Laws (Complete ✅)
- Chernoff upper/lower/two-sided bounds
- Hoeffding inequality for bounded random variables
- Adleman's Theorem: BPP ⊆ P/poly
- Sipser-Gács-Lautemann: BPP ⊆ Σ₂^p ∩ Π₂^p
- Schwartz-Zippel Lemma for polynomial identity testing
- Freivalds' Theorem for matrix multiplication verification
- Azuma-Hoeffding for martingale difference sequences
- Shannon's counting lower bound for circuits

## L5: Algorithms/Methods (Complete ✅)
- Miller-Rabin primality testing (ZPP)
- Freivalds' matrix verification (co-RP)
- Polynomial Identity Testing via Schwartz-Zippel (BPP)
- RP perfect matching via determinant
- Karger's randomized min-cut
- Papadimitriou's randomized 2-SAT walk
- Schöning's randomized 3-SAT walk
- Method of conditional probabilities (MAX-CUT, set balancing)
- Pairwise/k-wise independent hash construction

## L6: Canonical Problems (Complete ✅)
- Probabilistic Circuit-SAT
- PARITY circuit construction and evaluation
- MAJORITY circuit construction
- Approximate majority (Chernoff-based)
- Graph Non-Isomorphism (AM protocol)
- Global minimum cut

## L7: Applications (Complete ✅ — 3+ apps)
- Primality testing in cryptographic key generation
- Randomized quicksort expected-time analysis
- PCP verifier for proof checking
- Johnson-Lindenstrauss dimension reduction
- Matrix multiplication verification (distributed computing)

## L8: Advanced Topics (Partial ⚠️)
- Hardness vs randomness paradigm (NW generator concept)
- Impagliazzo's five worlds (documented)
- Circuit lower bounds for derandomization
- PCP Theorem connection

## L9: Research Frontiers (Partial — documented)
- Meta-complexity and derandomization
- Quantum vs classical randomness
- Fine-grained derandomization
- Natural proofs barrier and circuit lower bounds
