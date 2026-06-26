# Course Tree — mini-nisan-wigderson-prg

## Prerequisite Knowledge Dependencies

```
mini-nisan-wigderson-prg
│
├── mini-cook-levin-theorem (L4: NP-completeness framework)
│   └── Circuit SAT, reduction concepts
│
├── mini-p-np-np-completeness (L1-L2: Basic complexity classes)
│   └── P, NP, BPP definitions, circuit model
│
├── mini-boolean-circuits-model (L3: Circuit DAG model)
│   └── Gate types, evaluation, size/depth metrics
│
├── mini-time-hierarchy-theorem (L4: Hierarchy theorems)
│   └── EXP ≠ P, time constructible functions
│
├── mini-reductions-completeness (L2: Reduction types)
│   └── Turing reductions for hardness amplification
│
├── mini-hastad-lower-bounds (L8: Switching lemma)
│   └── AC0 lower bounds, PARITY ∉ AC0
│
└── mini-natural-proofs-barrier (L8: RR barrier)
    └── Natural properties, constructive arguments
```

## Internal Dependency Graph

```
nw_core.h  (defines: NWDesign, HardFunction, NWPRG, TruthTable)
    │
    ├── nw_designs.h  (defines: SetSystem, Expander; uses: NWDesign)
    │   └── nw_design.c  (RS design, random design, expander design)
    │
    ├── nw_hardness.h  (defines: HardnessAssumption, WorstToAverage; uses: TruthTable)
    │   └── nw_hardness.c  (Yao XOR, truth tables, Fourier, natural proofs)
    │
    ├── nw_circuits.h  (defines: Circuit; uses: NWDesign, TruthTable)
    │   └── nw_circuits.c  (PARITY, MAJORITY, switching lemma, GL predicate)
    │
    └── nw_prg.h  (includes all above; defines: PRGSecurity, StatisticalTestResult)
        └── nw_prg.c  (PRG construction, evaluation, hybrid arg, BPP sim, KDF)
```

## Knowledge Flow

```
Hardness Assumption (L1)
       │
       ├── Yao XOR Lemma (L4) ──> Hardness Amplification (L5)
       │
       ├── Worst-to-Average (L4) ──> Hard-Core Set (L5)
       │
       └── Combinatorial Design (L3)
              │
              ├── Reed-Solomon Construction (L5)
              ├── Random Design (L5)
              └── Expander Design (L8)
                     │
                     └── NW PRG Construction (L4)
                            │
                            ├── PRG Evaluation (L5)
                            ├── Hybrid Argument (L5)
                            ├── Distinguisher Analysis (L5)
                            │
                            └── Applications (L7)
                                   ├── BPP Derandomization
                                   ├── Key Derivation (KDF)
                                   ├── Statistical Testing
                                   └── IW97 (BPP = P)
```