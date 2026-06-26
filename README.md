# Mini Randomized & Derandomization

A collection of **from-scratch, zero-dependency C implementations** of university-level randomized computation and derandomization theory. Each module maps to MIT (and other top-tier university) courses, translating randomness-extraction, pseudorandom generation, and hardness-vs-randomness theorems into runnable C code.

## Sub-Modules

| Sub-Module | Topics | Key Courses |
|-----------|--------|-------------|
| [mini-bpp-rp-zpp-classes](mini-bpp-rp-zpp-classes/) | Randomized complexity classes (BPP/RP/ZPP/co-RP), probability amplification, Chernoff bounds, Schwartz-Zippel PIT, Freivalds verification | MIT 6.045, MIT 6.841, Stanford CS254 |
| [mini-derandomizing-space](mini-derandomizing-space/) | Space-bounded derandomization, Savitch's theorem, Immerman-Szelepcsényi (NL=coNL), Reingold's USTCON, Nisan's space PRG, branching programs | MIT 6.841/18.404, Stanford CS254, Berkeley CS278 |
| [mini-expander-graphs-construction](mini-expander-graphs-construction/) | Spectral expanders, Cheeger inequality, zigzag product, Ramanujan graphs (LPS), Margulis/Gabber-Galil constructions, expander mixing lemma | MIT 6.841, Stanford CS254, Berkeley CS278 |
| [mini-hardness-vs-randomness](mini-hardness-vs-randomness/) | Circuit lower bounds (Shannon, Håstad, Razborov), hardness amplification, worst-case to average-case reduction, Impagliazzo-Wigderson theorem | MIT 6.841, Stanford CS254, Berkeley CS278 |
| [mini-list-decodable-codes](mini-list-decodable-codes/) | List-decodable codes, Johnson/Griesmer bounds, Goldreich-Levin theorem, Reed-Solomon/Reed-Muller codes, folded Reed-Solomon | MIT 6.841, Stanford CS254, CMU 15-855 |
| [mini-nisan-wigderson-prg](mini-nisan-wigderson-prg/) | Nisan-Wigderson PRG construction, combinatorial designs, hybrid argument, Yao's XOR lemma, hardness amplification, BPP derandomization | MIT 6.841, Stanford CS254, Berkeley CS278 |
| [mini-pseudorandom-generators](mini-pseudorandom-generators/) | PRG constructions, next-bit unpredictability, NIST statistical tests, hardcore predicates, PRG composition, Yao's theorem | MIT 6.841, Stanford CS254, CMU 15-855, Princeton COS522 |
| [mini-randomness-extractors](mini-randomness-extractors/) | Randomness extractors, min-entropy, Leftover Hash Lemma, Trevisan extractor, dispersers, two-source extractors, entropy measures | MIT 6.841, Stanford CS254, Berkeley CS278, CMU 15-855 |

## Design Philosophy

- **Zero external dependencies** — pure C (C99/C11), only `libc` and `libm`
- **Self-contained modules** — each directory has its own `Makefile`, `include/`, `src/`, `examples/`, `demos/`, `docs/`
- **Theory-to-code mapping** — every module includes `docs/` with course-alignment notes and knowledge graphs
- **Derandomization pipeline** — each module connects to the broader hardness-vs-randomness framework (Arora-Barak Ch.19–20)

## Building

Each module is standalone. Navigate to a module directory and run:

```bash
cd mini-bpp-rp-zpp-classes
make all    # build everything
make test   # run tests
```

Requires **GCC** and **GNU Make**.

## Project Structure

```
mini-randomized-derandomization/
├── mini-bpp-rp-zpp-classes/       # Randomized complexity classes: BPP, RP, ZPP, co-RP
├── mini-derandomizing-space/      # Space-bounded derandomization: Savitch, NL=coNL, Reingold
├── mini-expander-graphs-construction/  # Spectral expanders, zigzag product, Ramanujan graphs
├── mini-hardness-vs-randomness/   # Circuit lower bounds, worst-to-average, IW theorem
├── mini-list-decodable-codes/     # List-decoding: Johnson bound, RS/RM codes, Goldreich-Levin
├── mini-nisan-wigderson-prg/      # NW PRG: combinatorial designs, hybrid argument, XOR lemma
├── mini-pseudorandom-generators/  # PRG foundations: NIST tests, hardcore predicates, composition
└── mini-randomness-extractors/    # Extractors: min-entropy, LHL, Trevisan, dispersers
```

## License

MIT
