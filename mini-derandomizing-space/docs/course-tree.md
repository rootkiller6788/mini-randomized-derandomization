# Course Tree ? mini-derandomizing-space

## Prerequisites
- Basic complexity theory (P, NP, polynomial time)
- Turing machines (deterministic and nondeterministic)
- Basic graph theory (connectivity, BFS/DFS)
- Basic probability theory (Chernoff bounds)

## This Module Provides
- Space complexity classes (L, NL, RL, BPL, PSPACE)
- Savitch theorem: NSPACE ? DSPACE(s^2)
- Immerman-Szelepcsenyi: NL = coNL
- Reingold theorem: USTCONN ? L
- Nisan PRG: BPL derandomization
- Expander graphs and zig-zag product

## Dependencies (modules this depends on)
- 0.mini-complexity-foundations/mini-p-np-np-completeness (basic classes)
- 0.mini-complexity-foundations/mini-time-hierarchy-theorem (diagonalization)

## Dependents (modules that depend on this)
- 2.mini-algorithmic-complexity/mini-pcp-theorem (derandomization techniques)
- 2.mini-algorithmic-complexity/mini-ip-pspace (interactive proofs)
- 2.mini-algorithmic-complexity/mini-approximation-algorithms (hardness of approximation)
- 1.mini-circuit-complexity/mini-natural-proofs-barrier (derandomization barriers)

## Learning Path
1. Start with space_derand.h/c for basic definitions
2. Study savitch.c for the recursive simulation technique
3. Study immerman_szelepcsenyi.c for inductive counting
4. Study nisan_prg.c for PRG construction
5. Study reingold_zigzag.c for expander-based derandomization
6. Read examples/ for end-to-end demonstrations
7. Read src/derandomizing_space.lean for formal verification
