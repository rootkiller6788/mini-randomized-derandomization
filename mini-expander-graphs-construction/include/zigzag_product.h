#ifndef ZIGZAG_H
#define ZIGZAG_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include "expander_core.h"

/* ------------------------------------------------------------------
 * Data Types
 * ------------------------------------------------------------------ */

/** ZigZagResult: holds input graphs G, H and the resulting product graph.
 *  G and H are borrowed (not owned); result is owned by this struct. */
typedef struct { ExpanderGraph *G; ExpanderGraph *H; ExpanderGraph *result; } ZigZagResult;

/* ------------------------------------------------------------------
 * Core Zig-Zag Operations
 * ------------------------------------------------------------------ */

/** Compute zig-zag product G (z) H. n_H must equal d_G. */
ZigZagResult zz_compute(const ExpanderGraph *G, const ExpanderGraph *H);

/** Rotation map: given (v,i), return (w,j) where edge i of v
 *  leads to w and this is edge j of w. */
size_t exp_rotation_map(const ExpanderGraph *g, size_t v, size_t i,
                         size_t *j);

/* ------------------------------------------------------------------
 * Spectral & Structural Properties
 * ------------------------------------------------------------------ */

/** Check if zig-zag product preserves expansion (RVW 2002 bound). */
bool zz_preserves_expansion(const ZigZagResult *z);

/** Degree reduction ratio: d_H^2 / d_G. */
double zz_degree_reduction(size_t d1, size_t d2);

/** New vertex count: n_G * d_G = n_G * n_H. */
size_t zz_new_size(size_t n1, size_t d2);

/** Theoretical upper bound on lambda_2(G(z)H). */
double zz_spectral_bound(const ExpanderGraph *G, const ExpanderGraph *H);

/** Empirical lambda_2 ratio (actual/theoretical). */
double zz_empirical_lambda(const ZigZagResult *z);

/* ------------------------------------------------------------------
 * Construction & Verification
 * ------------------------------------------------------------------ */

/** Build an explicit expander using one iteration of zig-zag + powering. */
bool zz_explicit_construction(size_t n, size_t d, ExpanderGraph **result);

/** Verify zig-zag product correctness (regularity, connectivity, size). */
bool zz_verify_construction(const ZigZagResult *z);

/* ------------------------------------------------------------------
 * Memory Management & Display
 * ------------------------------------------------------------------ */

/** Free zig-zag result (frees result graph, keeps G and H). */
void zz_free_result(ZigZagResult *z);

/** Print information about a zig-zag product to a stream. */
void zz_print_info(const ZigZagResult *z, FILE *fp);

#endif /* ZIGZAG_H */
