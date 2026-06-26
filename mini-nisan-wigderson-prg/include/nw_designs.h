/**
 * @file nw_designs.h
 * @brief Combinatorial design constructions for the NW PRG.
 *
 * A (k, m, l)-design is a family of m subsets of a universe of size k,
 * each of size exactly l, such that any two distinct subsets intersect in
 * at most log(m) elements.
 *
 * These designs are the combinatorial core of the Nisan-Wigderson
 * construction. The small intersection property ensures that the PRG's
 * output bits are "nearly independent."
 *
 * Construction methods:
 *   1. Reed-Solomon codes → designs (explicit, deterministic)
 *   2. Probabilistic method → random designs (existence proof)
 *   3. Expander graphs → designs (spectral approach)
 *
 * Reference:
 *   Nisan & Wigderson, J. Comput. Syst. Sci. 49(2), 1994
 *   Arora & Barak, Ch 20.2
 *   Jukna, "Boolean Function Complexity", Ch 11
 */

#ifndef NW_DESIGNS_H
#define NW_DESIGNS_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include "nw_core.h"

/* ──────────────────────────────────────────────
 * L1: Set System Definition
 * ────────────────────────────────────────────── */

/**
 * @brief A set system: collection of subsets of {0,...,n-1}.
 *
 * L1: Definition — Set system
 * L3: Mathematical Structure — Hypergraph
 */
typedef struct {
    size_t universe;      /**< n = |U| */
    size_t num_sets;      /**< m: number of subsets */
    size_t set_size;      /**< l: each set has size l (uniform) */
    size_t *elements;     /**< Flattened: sets[i*l + j] = j-th element of set i */
} SetSystem;

/**
 * @brief Expander graph for design construction.
 *
 * An (n, d, λ)-expander is a d-regular graph on n vertices
 * with second eigenvalue ≤ λ.
 *
 * L3: Mathematical Structure — Expander graph
 * L8: Advanced Topic — Spectral graph theory for derandomization
 */
typedef struct {
    size_t n;             /**< Number of vertices */
    size_t d;             /**< Degree */
    size_t *adj;          /**< Adjacency: adj[v*d + j] = j-th neighbor of v */
    double lambda;        /**< Second eigenvalue bound */
} Expander;

/* ──────────────────────────────────────────────
 * L5: Reed-Solomon Design Construction
 * ────────────────────────────────────────────── */

/**
 * @brief Construct a design using Reed-Solomon codes.
 *
 * Let F = GF(q) be a finite field. Universe = F × F (size q²).
 * For each polynomial p of degree ≤ d, define:
 *   S_p = { (x, p(x)) : x ∈ F }
 * This yields m = q^{d+1} sets of size q, with pairwise
 * intersection at most d.
 *
 * Choose parameters so that q = l and d ≤ intersection_bound.
 *
 * Theorem: For any ε > 0, there exist (k, m, l)-designs with
 * k = l², m = 2^{εl}, and intersection ≤ εl.
 *
 * L5: Algorithm — Reed-Solomon based design
 * L3: Mathematical Structure — Finite field GF(q)
 *
 * @param q               Field size (must be prime power ≤ 256)
 * @param degree          Maximum polynomial degree d
 * @return Constructed design, or NULL on error
 */
NWDesign *nw_reed_solomon_design(size_t q, size_t degree);

/**
 * @brief Construct a design using Reed-Solomon codes with explicit parameters.
 *
 * Given desired (k, m, l, bound), finds smallest prime power q
 * and degree d that satisfy the constraints.
 *
 * @param k   Universe size
 * @param m   Number of sets
 * @param l   Set size
 * @param bound  Intersection bound
 * @return Constructed design, or NULL if infeasible
 */
NWDesign *nw_reed_solomon_design_param(
    size_t k, size_t m, size_t l, size_t bound
);

/* ──────────────────────────────────────────────
 * L5: Random Design Construction (Probabilistic Method)
 * ────────────────────────────────────────────── */

/**
 * @brief Sample a random design using the probabilistic method.
 *
 * Each of the m sets is chosen independently and uniformly
 * from all l-element subsets of [k].
 *
 * Theorem: With high probability, the random design has
 * intersection ≤ log m when l = O(log m) and k is large enough.
 *
 * L5: Algorithm — Random design via probabilistic method
 * L4: Fundamental Law — Probabilistic method existence proof
 *
 * @param k   Universe size
 * @param m   Number of sets
 * @param l   Set size
 * @return Randomly generated design
 */
NWDesign *nw_random_design(size_t k, size_t m, size_t l);

/**
 * @brief Compute success probability that a random design
 * satisfies the small intersection property.
 *
 * Uses the union bound over all pairs of sets.
 *
 * P[|S_i ∩ S_j| > t] ≤ C(l, t+1) · (l/k)^{t+1} · (C(k-t-1, l-t-1) / C(k, l))
 *
 * L4: Fundamental Law — Probabilistic analysis of designs
 *
 * @param k      Universe size
 * @param m      Number of sets
 * @param l      Set size
 * @param bound  Intersection bound
 * @return Success probability ∈ [0, 1]
 */
double nw_random_design_success_prob(
    size_t k, size_t m, size_t l, size_t bound
);

/**
 * @brief Determine number of trials needed to find a good design.
 *
 * Uses the success probability to compute:
 * trials = ceil(log(1-confidence) / log(1 - p_success))
 *
 * @param k           Universe size
 * @param m           Number of sets
 * @param l           Set size
 * @param bound       Intersection bound
 * @param confidence  Desired confidence level
 * @return Number of random designs to try
 */
size_t nw_design_trials_needed(
    size_t k, size_t m, size_t l, size_t bound, double confidence
);

/* ──────────────────────────────────────────────
 * L8: Expander-based Design Construction
 * ────────────────────────────────────────────── */

/**
 * @brief Create an expander graph from known constructions.
 *
 * Uses the Margulis/Gabber-Galil explicit expander construction
 * for prime n.
 *
 * L8: Advanced Topic — Explicit expander constructions
 *
 * @param n  Number of vertices
 * @param d  Desired degree (must be achievable)
 * @return Constructed expander, or NULL
 */
Expander *nw_expander_create_margulis(size_t n, size_t d);

/**
 * @brief Create a design from a random walk on an expander.
 *
 * Takes length-l random walks from each vertex as sets.
 *
 * Theorem (Impagliazzo-Wigderson 1997): Expander walks yield
 * designs suitable for derandomization.
 *
 * L8: Advanced Topic — Expander-based derandomization
 *
 * @param e  Expander graph
 * @param l  Walk length (set size)
 * @return Constructed design
 */
NWDesign *nw_design_from_expander(const Expander *e, size_t l);

/**
 * @brief Compute the second eigenvalue bound of a d-regular expander.
 *
 * Alon-Boppana bound: λ ≥ 2√(d-1) - o(1) for any infinite family.
 * Ramanujan graphs achieve λ ≤ 2√(d-1).
 *
 * L8: Advanced Topic — Spectral expander bounds
 *
 * @param d  Degree
 * @return Lower bound on λ
 */
double nw_expander_alon_boppana(double d);

/**
 * @brief Free an expander.
 * @param e  Expander to free
 */
void nw_expander_free(Expander *e);

/* ──────────────────────────────────────────────
 * Set System Operations
 * ────────────────────────────────────────────── */

/**
 * @brief Create a uniform set system.
 *
 * @param n  Universe size
 * @param m  Number of sets
 * @param l  Uniform set size
 * @return Newly allocated SetSystem
 */
SetSystem *nw_setsystem_create(size_t n, size_t m, size_t l);

/**
 * @brief Get the j-th element of the i-th set.
 * @param ss  Set system
 * @param i   Set index
 * @param j   Element index within set
 * @return The element
 */
size_t nw_setsystem_get(const SetSystem *ss, size_t i, size_t j);

/**
 * @brief Set the j-th element of the i-th set.
 * @param ss   Set system
 * @param i    Set index
 * @param j    Element index
 * @param val  Value to set
 */
void nw_setsystem_set(SetSystem *ss, size_t i, size_t j, size_t val);

/**
 * @brief Verify small intersection property.
 * @param ss    Set system
 * @param bound Maximum allowed intersection
 * @return true if all pairs intersect in ≤ bound elements
 */
bool nw_setsystem_verify(const SetSystem *ss, size_t bound);

/**
 * @brief Compute intersection size of two sets.
 * @param ss  Set system
 * @param i   First set index
 * @param j   Second set index
 * @return |S_i ∩ S_j|
 */
size_t nw_setsystem_intersection(const SetSystem *ss, size_t i, size_t j);

/**
 * @brief Free a set system.
 * @param ss  Set system to free
 */
void nw_setsystem_free(SetSystem *ss);

/**
 * @brief Convert SetSystem to NWDesign.
 * @param ss    Source set system
 * @param bound Intersection bound for the design
 * @return Converted design
 */
NWDesign *nw_setsystem_to_design(const SetSystem *ss, size_t bound);

/**
 * @brief Convert NWDesign to SetSystem.
 * @param d  Source design
 * @return Converted set system
 */
SetSystem *nw_design_to_setsystem(const NWDesign *d);

/* ──────────────────────────────────────────────
 * Binary operations for computing binomial coefficients
 * ────────────────────────────────────────────── */

/**
 * @brief Compute binomial coefficient C(n, k).
 *
 * Uses multiplicative formula to avoid overflow where possible.
 *
 * @param n  Total items
 * @param k  Items to choose
 * @return C(n, k)
 */
double nw_binomial(size_t n, size_t k);

/**
 * @brief Compute log(C(n, k)) using Stirling's approximation.
 *
 * L3: Mathematical Structure — Asymptotic enumeration
 *
 * @param n  Total items
 * @param k  Items to choose
 * @return ln(C(n,k))
 */
double nw_log_binomial(size_t n, size_t k);

#endif /* NW_DESIGNS_H */
