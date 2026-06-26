# Course Tree ¡ª mini-hardness-vs-randomness

## Prerequisites
```
Boolean logic ¡ú Circuit complexity ¡ú Shannon lower bound
                                 ¡ú H?stad switching lemma
                                 ¡ú Razborov monotone bounds
                                 ¡ú Smolensky algebraic bounds
                                    ¡ý
                              Hardness definitions
                                    ¡ý
                              Yao's XOR Lemma ¡ú Worst-to-average
                              Impagliazzo Hardcore Lemma
                                    ¡ý
                              Nisan-Wigderson PRG construction
                                    ¡ý
                              Impagliazzo-Wigderson Theorem
                              (BPP = P if EXP has hard functions)
                                    ¡ý
                              Applications:
                              - Cryptography (OWF)
                              - Derandomization
                              - Complexity classification
```

## Internal Dependencies
```
circuit_lower_bounds.{h,c}
    ¡ý
hardness_randomness.{h,c}  ¡û  worst_to_average.{h,c}
    ¡ý                              ¡ý
derandomization_via_hardness.{h,c}
    ¡ý
iw_theorem.c  ¡û  applications.c
    ¡ý
circuit_simulation.c
```

## External Dependencies (Cross-Module)
- mini-cook-levin-theorem: P, NP definitions needed for BPP context
- mini-bpp-rp-zpp-classes: BPP definition and probability amplification
- mini-nisan-wigderson-prg: NW PRG construction details
- mini-expander-graphs-construction: Expander walks for derandomization
- mini-pseudorandom-generators: PRG security definitions
- mini-circuits-model (circuit complexity): Circuit class definitions

## Key Theorems in Dependency Order
1. Shannon (1949): Counting lower bound ¡ú existential hardness
2. Lupanov (1958): Upper bound ¡ú Shannon is tight
3. Yao (1982): XOR lemma ¡ú hardness amplification
4. H?stad (1986): Switching lemma ¡ú PARITY ? AC0
5. Razborov (1985): Monotone bounds ¡ú CLIQUE hardness
6. Smolensky (1987): Algebraic bounds ¡ú MODp ? AC0[q]
7. Goldreich-Levin (1989): Hardcore predicate ¡ú PRG building block
8. Nisan-Wigderson (1994): Hardness ¡ú PRG
9. Impagliazzo (1995): Five worlds + Hardcore lemma
10. Impagliazzo-Wigderson (1997): P = BPP under hardness