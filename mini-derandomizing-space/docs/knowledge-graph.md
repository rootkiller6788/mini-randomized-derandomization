# Knowledge Graph ? mini-derandomizing-space

## L1: Definitions
| Entry | Status | Location |
|-------|--------|----------|
| SpaceTM (space-bounded Turing machine) | Complete | include/space_derand.h |
| TMConfig (configuration) | Complete | include/space_derand.h |
| SpaceClass (L, NL, RL, BPL, PSPACE, etc.) | Complete | include/space_derand.h |
| ConfigGraph (configuration graph) | Complete | include/space_derand.h |
| SpaceConstructible (constructibility) | Complete | include/space_derand.h |
| PairwiseHashFamily (pairwise independent hash) | Complete | include/space_derand.h |
| SpacePRG (pseudorandom generator) | Complete | include/space_derand.h |
| RegularGraph (d-regular graph) | Complete | include/space_derand.h |
| SpectralData (eigenvalue data) | Complete | include/space_derand.h |
| BranchingProgram (layered DAG, L/poly) | Complete | include/space_derand.h |
| LogspaceFn (logspace-computable function) | Complete | include/logspace_derand.h |
| NisanHashDesc (Nisan hash descriptor) | Complete | include/nisan_prg.h |

## L2: Core Concepts
| Entry | Status | Location |
|-------|--------|----------|
| Class containments (L?NL?P?PSPACE) | Complete | src/space_derand.c sd_class_contains |
| Configuration graph reachability | Complete | src/space_derand.c |
| Space constructibility | Complete | src/space_derand.c sd_log_constructible |
| Error amplification (Chernoff bound) | Complete | src/space_derand.c sd_error_amplify |
| Logspace derandomization | Complete | src/logspace_derand.c |
| Deterministic vs nondeterministic space | Complete | src/savitch.c |
| Randomized space (RL, BPL) | Complete | src/nisan_prg.c |
| Simulation equivalence | Complete | src/space_derand.c sd_verify_simulation |

## L3: Mathematical Structures
| Entry | Status | Location |
|-------|--------|----------|
| Boolean matrix representation (A*x + b over GF(2)) | Complete | src/space_derand.c sd_hash_eval |
| Regular graph rotation maps | Complete | src/space_derand.c sd_regular_graph_create |
| Normalized adjacency matrix | Complete | src/reingold_zigzag.c rz_adjacency_matrix |
| Eigenvalue decomposition (power iteration) | Complete | src/space_derand.c sd_power_iteration |
| Branching program evaluation (DP) | Complete | src/space_derand.c sd_bp_accept_probability |
| Configuration space enumeration | Complete | src/savitch.c |
| Graph products (zig-zag, tensor, replacement) | Complete | src/space_derand.c, src/reingold_zigzag.c |

## L4: Fundamental Theorems
| Theorem | Status | Location |
|---------|--------|----------|
| Savitch: NSPACE(s) ? DSPACE(s^2) | Complete | src/savitch.c, src/derandomizing_space.lean |
| Immerman-Szelepcsenyi: NL = coNL | Complete | src/immerman_szelepcsenyi.c, src/derandomizing_space.lean |
| Reingold: USTCONN ? L (SL = L) | Complete | src/logspace_derand.c, src/derandomizing_space.lean |
| Nisan PRG: seed O(S log T) | Complete | src/nisan_prg.c, src/derandomizing_space.lean |
| Alon-Boppana bound | Complete | src/reingold_zigzag.c |

## L5: Algorithms
| Algorithm | Status | Location |
|-----------|--------|----------|
| Savitch recursive simulation | Complete | src/savitch.c |
| Immerman-Szelepcsenyi inductive counting | Complete | src/immerman_szelepcsenyi.c |
| Reingold expander-based USTCONN | Complete | src/logspace_derand.c, src/reingold_zigzag.c |
| Nisan PRG construction and evaluation | Complete | src/nisan_prg.c |
| Zig-zag product construction | Complete | src/space_derand.c |
| Graph powering for expansion boost | Complete | src/space_derand.c |
| Spectral gap computation (power iteration) | Complete | src/space_derand.c |
| Pairwise independent hash family | Complete | src/space_derand.c |
| BFS/DFS graph traversal | Complete | src/space_derand.c |
| Union-find component counting | Complete | src/logspace_derand.c |

## L6: Canonical Problems
| Problem | Status | Location |
|---------|--------|----------|
| USTCONN (undirected s-t connectivity) | Complete | src/logspace_derand.c, examples/example1.c |
| STCONN (directed s-t connectivity) | Complete | src/logspace_derand.c |
| BPL-acceptance | Complete | src/nisan_prg.c |
| RL-acceptance | Complete | src/space_derand.c |
| Connectivity (graph) | Complete | src/space_derand.c |
| Component counting | Complete | src/logspace_derand.c |
| Shortest path (unweighted) | Complete | src/space_derand.c |

## L7: Applications
| Application | Status | Location |
|-------------|--------|----------|
| Derandomized network connectivity (sensor nets) | Complete | src/logspace_derand.c sd_network_connectivity_check |
| Model checking safety properties | Complete | src/logspace_derand.c sd_safety_property_check |
| Expander-based routing networks | Complete | src/reingold_zigzag.c rz_routing_network_build |
| BPL derandomization for cryptographic protocols | Complete | src/nisan_prg.c nisan_derandomize_bpl |
| Error-correcting codes from expanders | Complete | src/reingold_zigzag.c rz_expander_code_check |

## L8: Advanced Topics
| Topic | Status | Location |
|-------|--------|----------|
| Expander mixing lemma | Complete | src/reingold_zigzag.c rz_mixing_lemma_check |
| Ramanujan graphs | Complete | src/reingold_zigzag.c rz_is_ramanujan |
| Derandomized graph squaring | Complete | src/logspace_derand.c sd_derandomized_square |
| PRG stretch analysis | Complete | src/nisan_prg.c nisan_stretch_factor |
| Statistical distance computation | Complete | src/nisan_prg.c nisan_statistical_distance |

## L9: Research Frontiers
| Topic | Status | Location |
|-------|--------|----------|
| High-dimensional expanders | Partial | src/reingold_zigzag.c rz_high_dimensional_expander_test |
| RL vs L (open problem) | Partial | src/derandomizing_space.lean RL_equals_L_open |
| BPL derandomization frontier | Partial | src/derandomizing_space.lean |
