# Course Alignment — mini-expander-graphs-construction

## Nine-School Curriculum Mapping

| School | Course | Expander Topics Covered | Implementation |
|--------|--------|------------------------|---------------|
| **MIT** | 6.841 Advanced Complexity | Zig-zag product, SL=L, derandomization via expanders | `zz_compute()`, `exp_sl_connectivity()`, `exp_error_reduction()` |
| **Stanford** | CS254 Computational Complexity | Spectral expanders, Cheeger inequality, mixing time | `exp_cheeger_constant()`, `exp_random_walk_mix_time()` |
| **Berkeley** | CS278 Computational Complexity | LPS Ramanujan graphs, PGL(2,q) construction | `exp_ramanujan_lps()` |
| **CMU** | 15-855 Graduate Complexity | Expander mixing lemma, expander codes | `exp_mixing_lemma_bound()`, `exp_code_encode()` |
| **Princeton** | COS 522 Computational Complexity | Expanders for cryptography, extractors, PRGs | `exp_extractor_graph()`, `exp_prg_generate()` |
| **Caltech** | CS 151 Complexity Theory | Spectral graph theory, Fiedler vector, partitioning | `exp_fiedler_vector()`, `exp_spectral_partition()` |
| **Cambridge** | Part III Advanced Complexity | Cayley graphs, group-theoretic expander constructions | `exp_cayley_graph()`, `exp_margulis_gabber_galil()` |
| **Oxford** | Advanced Complexity Theory | Lossless expanders, high-dimensional expanders | `exp_is_lossless_expander()` (heuristic) |
| **ETH** | 263-4650 Advanced Complexity | Deterministic log-space algorithms, Reingold's theorem | `exp_sl_path_length()`, `zz_explicit_construction()` |

## Cross-Cutting Themes
- **Derandomization**: MIT, Stanford, Princeton, ETH — error reduction via expander walks
- **Spectral Methods**: Stanford, Caltech, Cambridge — eigenvalue-based analysis
- **Explicit Constructions**: Berkeley, Cambridge, ETH — Margulis, LPS, zig-zag
- **Applications**: CMU, Princeton, Oxford — codes, cryptography, extractors