# Knowledge Graph — mini-randomness-extractors

## L1: Definitions

| # | Definition | C Implementation | Lean |
|---|-----------|-----------------|------|
| 1 | Extractor (k,eps)-extractor Ext : {0,1}^n x {0,1}^d -> {0,1}^m | `Extractor` struct (extractor_core.h) | `Extractor` structure |
| 2 | Weak Random Source (k-source) | `WeakSource` struct | `KSource` structure |
| 3 | Seeded Extractor | `SeededExt` struct with function pointer | – |
| 4 | Two-Source Extractor | `TwoSourceExt` struct | – |
| 5 | Disperser | `Disperser` struct (disperser_core.h) | – |
| 6 | Oblivious Disperser | `ObliviousDisp` struct | – |
| 7 | Min-entropy H_inf | `ext_min_entropy` | `min_entropy_uniform` theorem |
| 8 | Statistical Distance Delta | `ext_statistical_distance` | `stat_distance_bounded` theorem |
| 9 | Universal Hash Family | `UniversalHash` struct | `carter_wegman_universal` theorem |
| 10 | Extractor Graph | `ExtractorGraph` struct | – |

## L2: Core Concepts

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Min-entropy = -log2(max probability) | `ext_min_entropy_from_freq` |
| 2 | Shannon entropy H = -sum p log p | `ext_shannon_entropy` |
| 3 | Renyi entropy H_alpha | `ext_renyi_entropy` |
| 4 | Total variation / statistical distance | `ext_statistical_distance` |
| 5 | Hellinger distance | `ext_hellinger_distance` |
| 6 | KL divergence D_KL | `ext_kl_divergence` |
| 7 | Collision probability CP(X) | `collisionProbability` (Lean) |
| 8 | Entropy deficit (gap to uniform) | `ext_entropy_deficit` |
| 9 | Data processing inequality | See stat_distance_bounded |
| 10 | Pinsker inequality | Example 1 demonstration |

## L3: Mathematical Structures

| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | Probability simplex over {0,1}^n | Frequency table operations |
| 2 | GF(2) field operations | `ext_gf2_add`, `ext_gf2_mul` |
| 3 | Bipartite expander graph | `ExtractorGraph` + spectral bounds |
| 4 | Carter-Wegman hash family over prime fields | `uh_create_carter_wegman` |
| 5 | Pairwise independent hash over GF(2) | `uh_create_pairwise_independent` |
| 6 | Combinatorial designs | `design_create` |
| 7 | Reed-Solomon codes (for Ta-Shma) | Section 6 of extractor_constructions.c |
| 8 | Bit-vector algebra | `ext_bits_pack/unpack/popcount` |

## L4: Fundamental Laws

| # | Theorem | Statement | Verification |
|---|---------|-----------|-------------|
| 1 | Leftover Hash Lemma (ILL 1989) | If H_inf(X) >= k, Delta((H,H(X)),(H,U_m)) <= 2^{-(k-m)/2} | `ext_error_bound`, `leftover_hash_lemma` (Lean) |
| 2 | RTD Bound (2000) | m <= k + 2*log(eps) + 2 | `ext_output_len` |
| 3 | NZ Seed Bound (1996) | d >= log(n-k) + 2*log(1/eps) - O(1) | `ext_seed_len` |
| 4 | Chernoff Inequality | Pr[|p_hat - p| > eps] <= 2*exp(-2*n*eps^2) | `ext_verify_min_entropy` |
| 5 | H_inf <= H (Shannon) | Min-entropy never exceeds Shannon | Test + `min_entropy_le_shannon` |
| 6 | Pinsker Inequality | Delta(P,Q) <= sqrt(D_KL(P||Q) / 2) | Example verification |
| 7 | Data Processing | Delta(f(P),f(Q)) <= Delta(P,Q) | Implicit in structure |
| 8 | Gibbs Inequality | D_KL(P||Q) >= 0 | Test `ext_kl_divergence` |
| 9 | Triangle inequality for Delta | Delta(P,R) <= Delta(P,Q) + Delta(Q,R) | Test verified |
| 10 | Disperser from extractor | Every extractor is a disperser | `disp_from_extractor` |

## L5: Algorithms/Methods

| # | Algorithm | Implementation |
|---|-----------|---------------|
| 1 | Carter-Wegman universal hashing | `uh_create_carter_wegman`, `uh_eval` |
| 2 | Leftover Hash Lemma extractor | `ext_leftover_hash` |
| 3 | Trevisan extractor (NW + designs) | `ext_trevisan` |
| 4 | Raz block-and-condense extractor | `ext_raz` |
| 5 | Ta-Shma code-based extractor | `ext_ta_shma` |
| 6 | Chor-Goldreich inner-product 2-source | `ext_two_source_chor_goldreich` |
| 7 | Bourgain sum-product 2-source | `ext_two_source_bourgain` |
| 8 | NW PRG from hard predicate | `ext_nw_generator` |
| 9 | Block-source extractor | `ext_block_source` |
| 10 | Extractor graph construction | `eg_construct` |
| 11 | Disperser via polynomial evaluation | `disp_eval` |
| 12 | Seeded extractor verification (chi-squared) | `ext_seeded_verify` |

## L6: Canonical Problems

| # | Problem | Example/Demo |
|---|---------|-------------|
| 1 | Extract from biased coin source | Example 1, Part 2 |
| 2 | Estimate statistical distance to uniform | Example 1, Part 4 |
| 3 | Verify min-entropy of a source | Tests, `ext_verify_min_entropy` |
| 4 | Build extractor for given parameters | Example 2 |
| 5 | 2-source extraction from independent sources | Example 2, Part 4 |
| 6 | Derandomize BPP algorithm | Example 3, Part 1 |
| 7 | Key derivation from DH output | Example 3, Part 2 |
| 8 | Privacy amplification after QKD | Example 3, Part 3 |

## L7: Applications

| # | Application | Implementation | Domain |
|---|------------|---------------|--------|
| 1 | BPP derandomization (NW/IW) | `ext_derandomize_bpp` | Complexity theory |
| 2 | NW hardness-vs-randomness | `ext_bpp_nw_derandomize` | Derandomization |
| 3 | Cryptographic key derivation | `ext_key_derivation` | Cryptography |
| 4 | Privacy amplification (BBR88) | `ext_privacy_amplification` | QKD / Info theory |
| 5 | Explicit extractor simulation | `ext_simulate_explicit_extractor` | Empirical validation |
| 6 | Statistical test battery | `ext_statistical_test_battery` | NIST-style testing |

## L8: Advanced Topics

| # | Topic | Implementation |
|---|-------|---------------|
| 1 | Trevisan extractor (optimal seed) | `ext_trevisan` |
| 2 | Bourgain 2-source (sum-product) | `ext_two_source_bourgain` |
| 3 | Extractor graphs & expansion | `eg_construct`, spectral bounds |
| 4 | Network extractor protocols (KRV09) | `ext_network_extractor_simulate` |
| 5 | Multi-source fault tolerance | `ext_multi_source_fault_tolerant` |
| 6 | Quantum-proof extraction | `ext_quantum_proof_bound` |

## L9: Research Frontiers

| # | Topic | Documentation |
|---|-------|--------------|
| 1 | Quantum-proof extractors | Security bound analysis |
| 2 | Network extractors (KRV 2009) | Protocol simulation |
| 3 | Multi-source extractors | Fault-tolerant implementation |
| 4 | Non-malleable extractors | Documented in course-tree |
| 5 | Extractors for space-bounded computation | NZ paradigm referenced |
