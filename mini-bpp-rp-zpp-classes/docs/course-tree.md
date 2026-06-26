# Course Tree — Prerequisites for BPP/RP/ZPP Classes

## Direct Prerequisites
- **mini-p-np-np-completeness** — P vs NP context for randomized classes
- **mini-time-hierarchy-theorem** — Time hierarchy for EXP containment
- **mini-boolean-circuits-model** — Circuit model for BPP ⊆ P/poly

## Conceptual Prerequisites
- Probability theory (basic): random variables, expectation, variance
- Turing machine model: deterministic, non-deterministic
- Asymptotic notation: O, Ω, Θ
- Basic combinatorics: binomial coefficients, union bounds

## Dependencies
```
mini-p-np-np-completeness
    └── mini-bpp-rp-zpp-classes (this module)
            ├── mini-pseudorandom-generators
            │       └── mini-nisan-wigderson-prg
            ├── mini-randomness-extractors
            ├── mini-hardness-vs-randomness
            ├── mini-derandomizing-space
            ├── mini-expander-graphs-construction
            └── mini-list-decodable-codes
```

## Learning Path
1. Start with deterministic P and NP definitions
2. Understand why randomness helps (examples: identity testing, perfect matching)
3. Define BPP, RP, ZPP formally
4. Learn Chernoff/Hoeffding bounds for error analysis
5. Study BPP = P conditional on hardness assumptions
6. Explore advanced topics: interactive proofs, PCP, extractors
