# Knowledge Graph — mini-nisan-wigderson-prg

## L1: Core Definitions

| # | Definition | Implementation |
|---|-----------|---------------|
| 1 | NW combinatorial (k,m,l)-design | `NWDesign` in nw_core.h, nw_design.c |
| 2 | Hard function (circuit-size hardness) | `HardFunction` in nw_core.h |
| 3 | Boolean circuit over {AND,OR,NOT} | `Circuit` in nw_circuits.h |
| 4 | Pseudorandom generator G: {0,1}^k->{0,1}^m | `NWPRG` in nw_core.h |
| 5 | Truth table (size 2^n) | `TruthTable` in nw_core.h, nw_hardness.c |
| 6 | Distinguisher / next-bit predictor | `Distinguisher`, `NextBitPredictor` in nw_core.h |
| 7 | Set system (uniform hypergraph) | `SetSystem` in nw_designs.h |
| 8 | Expander graph (n,d,lambda)-expander | `Expander` in nw_designs.h |
| 9 | Circuit lower bound | `CircuitBound` in nw_core.h |
| 10 | Hardness assumption | `HardnessAssumption` in nw_hardness.h |

## L2: Core Concepts

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Hardness vs randomness tradeoff | nw_prg.c (stretch computation) |
| 2 | Derandomization via hard functions | nw_derandomize.c, nw_prg.c |
| 3 | Computational indistinguishability | nw_prg.c (distinguisher advantage) |
| 4 | Next-bit unpredictability = PRG | nw_prg.c (hybrid argument) |
| 5 | Stretch parameter (m/k ratio) | nw_seed_length, nw_output_length |
| 6 | Hardness amplification | nw_hardness.c (Yao XOR Lemma) |
| 7 | PRG construction strategies | NWPRGMode in nw_prg.h |
| 8 | BPP / derandomization | nw_bpp_derandomize in nw_prg.c |

## L3: Mathematical Structures

| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | Set system with intersection constraints | nw_design.c |
| 2 | Finite field GF(q) for Reed-Solomon | nw_design.c (_mod_add, _poly_eval) |
| 3 | Boolean circuit DAG | nw_circuits.c (gate types, topological order) |
| 4 | Truth table representation (2^n bits) | nw_hardness.c |
| 5 | Fourier analysis of boolean functions | nw_fourier_coefficient, nw_variable_influence |
| 6 | Expander graph (Margulis construction) | nw_design.c |
| 7 | Binomial coefficients (asymptotic) | nw_binomial, nw_log_binomial |
| 8 | Hypergeometric distribution | nw_random_design_success_prob |
| 9 | Statistical distance | nw_statistical_distance |
| 10 | Hamming weight / Boolean XOR | nw_hamming_weight, nw_bool_xor |

## L4: Fundamental Laws

| # | Theorem | Implementation |
|---|---------|---------------|
| 1 | **Nisan-Wigderson Theorem** (Hardness => PRG) | nw_prg.c, nw_theorem.lean |
| 2 | **Yao's XOR Lemma** (Hardness amplification) | nw_yao_xor_amplify, nw_yao_xor_lemma_full |
| 3 | **Shannon's Counting Lower Bound** (1949) | nw_shannon_lower_bound |
| 4 | **Lupanov's Upper Bound** (1958) | nw_lupanov_upper_bound |
| 5 | **Impagliazzo's Hard-Core Set Theorem** | nw_worst_to_average |
| 6 | **Impagliazzo-Wigderson Theorem** (BPP=P if E hard) | nw_iw_derandomize |
| 7 | **Hastad's Switching Lemma** (AC0 lower bounds) | nw_hastad_applies, nw_switching_lemma_bound |
| 8 | **Razborov-Rudich Natural Proofs Barrier** | nw_natural_proofs_barrier |
| 9 | **Yao's Distinguisher-to-Predictor Theorem** | nw_distinguisher_to_predictors |
| 10 | **Goldreich-Levin Theorem** (Hard-core predicate) | nw_goldreich_levin_hard_predicate |

## L5: Algorithms/Methods

| # | Algorithm | Implementation |
|---|-----------|---------------|
| 1 | Reed-Solomon design construction | nw_reed_solomon_design |
| 2 | Random design (probabilistic method) | nw_random_design |
| 3 | Expander-based design construction | nw_design_from_expander |
| 4 | NW PRG evaluation (G(x)_i = f(x|_Si)) | nw_prg_evaluate, nw_prg_eval_range |
| 5 | Optimal PRG parameter selection | nw_prg_construct_optimal |
| 6 | Hybrid argument for NW proof | nw_hybrid_position, nw_hybrid_argument |
| 7 | Distinguisher advantage estimation (MC) | nw_prg_distinguisher_advantage |
| 8 | BPP derandomization via NW PRG | nw_bpp_derandomize |
| 9 | Statistical testing suite (NIST-like) | nw_prg_statistical_tests |
| 10 | Key derivation via PRG | nw_prg_key_derivation |
| 11 | Randomized algorithm simulation | nw_prg_simulate_randomized |
| 12 | Circuit compilation from truth table | nw_compile_truth_table |
| 13 | Hard-core set construction | nw_worst_to_average |
| 14 | Hardness assumption validation | nw_hardness_assumption_valid |
| 15 | Truth table creation/evaluation | nw_truth_table_create, nw_truth_table_eval |

## L6: Canonical Problems

| # | Problem | Implementation |
|---|---------|---------------|
| 1 | **PARITY** function (PARITY ∈ NC1, ∉ AC0) | nw_circuit_parity |
| 2 | **MAJORITY** function (MAJ ∈ TC0, ∉ AC0) | nw_circuit_majority |
| 3 | **MOD-p** function (Razborov-Smolensky) | nw_circuit_mod_p |
| 4 | Inner Product (Goldreich-Levin) | nw_goldreich_levin_hard_predicate |
| 5 | PRG Distinguishing Problem | nw_prg_distinguisher_advantage |
| 6 | Circuit Lower Bound (counting argument) | nw_shannon_circuit_lower_bound |

## L7: Applications

| # | Application | Implementation |
|---|------------|---------------|
| 1 | BPP derandomization (SUBEXP/P) | nw_bpp_derandomize, nw_iw_derandomize |
| 2 | Cryptographic key derivation (KDF) | nw_prg_key_derivation |
| 3 | Randomized algorithm simulation | nw_prg_simulate_randomized |
| 4 | PRG parameter optimization for BPP | nw_prg_for_BPP |
| 5 | Statistical testing of PRG output | nw_prg_statistical_tests |

## L8: Advanced Topics

| # | Topic | Implementation |
|---|-------|---------------|
| 1 | Håstad's Switching Lemma | nw_switching_lemma_bound, nw_hastad_applies |
| 2 | Natural Proofs Barrier | nw_natural_proofs_barrier |
| 3 | Expander-based derandomization | nw_expander_create_margulis |
| 4 | Razborov-Rudich threshold | nw_razborov_rudich_threshold |
| 5 | AC0 PARITY lower bound | nw_ac0_parity_lower_bound |
| 6 | Spectral expander bounds | nw_expander_alon_boppana |
| 7 | Worst-case to average-case reduction | nw_worst_to_average |

## L9: Research Frontiers

| # | Topic | Status |
|---|-------|--------|
| 1 | Meta-complexity (MCSP, MKTP) | Documented |
| 2 | Derandomizing BPP unconditionally | Documented |
| 3 | Circuit lower bounds beyond natural proofs | Documented |
| 4 | Quantum derandomization | Documented |
| 5 | Fine-grained derandomization | Documented |