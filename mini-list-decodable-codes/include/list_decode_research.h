/*
 * list_decode_research.h — L9 Research Frontiers in List-Decodable Codes
 *
 * Implements connections between list-decodable codes and
 * active research areas in computational complexity theory.
 *
 * L9 Topics:
 *   1. Approximate (relaxed) list-decoding
 *   2. Subspace evasive sets for capacity-achieving codes
 *   3. Quantum CSS codes from list-decodable classical codes
 *   4. Meta-complexity: MCSP and list-decoding
 *   5. Locally testable codes (LTCs) as list-decodable codes
 *
 * References:
 *   - Guruswami, "Relaxed Decoding and Approximate List Decoding" (2003)
 *   - Dvir & Lovett, "Subspace Evasive Sets" (STOC 2012)
 *   - Calderbank, Rains, Shor, Sloane, "Quantum Error Correction
 *     via Codes over GF(4)" (IEEE Trans. IT 1998)
 *   - Guruswami & Kopparty, "List-Decoding of LTCs" (unpublished note)
 *   - Kabanets & Cai, "Circuit Minimization Problem" (STOC 2000)
 */
#ifndef LD_RESEARCH_H
#define LD_RESEARCH_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

/* ------------------------------------------------------------------ */
/*  L9.1 — Approximate (Relaxed) List-Decoding                        */
/* ------------------------------------------------------------------ */

/**
 * Approximate list-decoding parameters.
 *
 * In relaxed list-decoding (Guruswami 2003), the decoder is allowed
 * to output codewords that are "approximately" in the code — i.e.,
 * within a relaxed Hamming ball of radius (1-α)·d of some codeword,
 * where d is the code's minimum distance.
 *
 * This relaxation allows decoding up to a larger error radius.
 */
typedef struct {
    size_t n;            /* block length */
    size_t k;            /* message length */
    size_t q;            /* alphabet size */
    double alpha;        /* relaxation parameter: radius = (1-alpha)*d */
    double error_frac;   /* target error fraction for list-decoding */
    size_t max_list;     /* maximum output list size */
} RelaxedDecodeParams;

/**
 * Relaxed list-decoding: find all messages whose encoding is within
 * (1-alpha)*d Hamming distance of the received word.
 *
 * @param received   Received word (length n over [q]).
 * @param rdp        Relaxed decoding parameters.
 * @param decoded    Output: array of candidate messages (flat).
 * @return Number of candidates found.
 */
size_t ld_relaxed_list_decode(const size_t *received,
                               const RelaxedDecodeParams *rdp,
                               size_t **decoded);

/**
 * Compute the achievable relaxed decoding radius.
 *   δ_relaxed(n, k, d, α) = (1-α)·d / n
 */
double ld_relaxed_decoding_radius(size_t n, size_t k, size_t d, double alpha);

/**
 * Test if relaxed list-decoding is possible for given parameters.
 * Uses the relaxed Johnson bound.
 */
bool ld_relaxed_is_decodable(const RelaxedDecodeParams *rdp);

/* ------------------------------------------------------------------ */
/*  L9.2 — Subspace Evasive Sets                                      */
/* ------------------------------------------------------------------ */

/**
 * A subspace evasive set S ⊂ GF(q)^n has the property that
 * for any k-dimensional subspace H ⊂ GF(q)^n,
 *   |S ∩ H| ≤ ℓ.
 *
 * These sets are crucial for constructing capacity-achieving
 * list-decodable codes (Guruswami-Xing 2012, Dvir-Lovett 2012).
 *
 * Construction: S = { (x, f(x)) : x ∈ GF(q)^(n-k) }
 * where f is a polynomial of appropriate degree.
 */
typedef struct {
    size_t n;            /* ambient dimension */
    size_t k;            /* subspace dimension to evade */
    size_t ell;          /* maximum intersection size */
    size_t q;            /* field size */
    size_t size;         /* |S| */
    size_t *points;      /* flat array: size × n, points[i*n + j] */
} SubspaceEvasiveSet;

/**
 * Construct a subspace evasive set using the polynomial method.
 * f(x) = (x_1, ..., x_{n-k}) → GF(q)^k is chosen to make
 * the graph of f subspace-evasive.
 *
 * Complexity: O(q^{n-k}·n).
 */
SubspaceEvasiveSet *ses_construct_polynomial(size_t n, size_t k,
                                              size_t q, size_t ell);

/**
 * Verify the subspace evasive property: for a random k-dimensional
 * subspace H, check that |S ∩ H| ≤ ℓ.
 *
 * @param num_trials Number of random subspaces to test.
 * @return true if all tested subspaces satisfy the bound.
 */
bool ses_verify_evasive(const SubspaceEvasiveSet *ses,
                         size_t num_trials);

/**
 * Free a subspace evasive set.
 */
void ses_free(SubspaceEvasiveSet *ses);

/**
 * Compute the existential bound on subspace evasive set size.
 * Theorem (Dvir-Lovett): There exists an (n,k,ell,q)-SES of size
 * at least q^{n-k-ell+1}.
 */
size_t ses_existential_size_bound(size_t n, size_t k, size_t q, size_t ell);

/* ------------------------------------------------------------------ */
/*  L9.3 — Quantum CSS Codes from Classical List-Decodable Codes      */
/* ------------------------------------------------------------------ */

/**
 * A CSS (Calderbank-Shor-Steane) quantum code is constructed from
 * two classical linear codes C1 ⊆ C2 over GF(2) with parameters
 * [n, k1, d1] and [n, k2, d2] respectively.
 *
 * The resulting quantum code encodes k2-k1 logical qubits into
 * n physical qubits with distance min(d1, d2^⊥).
 *
 * If C1 and C2 are list-decodable, the quantum code inherits
 * list-decodability properties (Leung-Smith 2009).
 */
typedef struct {
    size_t n;            /* physical qubits (codeword length) */
    size_t k_logical;    /* logical qubits = k2 - k1 */
    size_t d_x;          /* X-distance (from C1) */
    size_t d_z;          /* Z-distance (from C2^⊥) */
    /* The two classical codes */
    size_t k1;           /* dimension of C1 */
    size_t k2;           /* dimension of C2 (C1 ⊂ C2) */
    size_t *gen_c1;      /* generator for C1 (k1 × n, flat) */
    size_t *gen_c2;      /* generator for C2 (k2 × n, flat) */
} QuantumCSSCode;

/**
 * Construct a quantum CSS code from two list-decodable classical codes.
 *
 * The classical codes C1 and C2 must satisfy C1 ⊂ C2.
 * Their list-decodability determines the quantum code's
 * resilience against adversarial noise.
 *
 * @param gen_c1  Generator matrix for C1 (k1 × n, flat).
 * @param gen_c2  Generator matrix for C2 (k2 × n, flat).
 * @return A quantum CSS code, or NULL if C1 ⊄ C2.
 */
QuantumCSSCode *qcss_create_from_classical(size_t n, size_t k1, size_t k2,
                                            const size_t *gen_c1,
                                            const size_t *gen_c2);

/**
 * Quantum list-decoding: given a syndrome, find all error patterns
 * of weight ≤ t that are consistent with the observed syndrome.
 *
 * This reduces to list-decoding the classical codes C1 and C2^⊥.
 *
 * @param syndrome_x  X-syndrome (from C1 parity checks).
 * @param syndrome_z  Z-syndrome (from C2^⊥ parity checks).
 * @param t           Maximum error weight.
 * @param errors_out  Output error candidates (L × n, flat).
 * @param max_list    Maximum list size.
 * @return Number of candidate error patterns.
 */
size_t qcss_list_decode_errors(const QuantumCSSCode *qcss,
                                const size_t *syndrome_x,
                                const size_t *syndrome_z,
                                size_t t,
                                size_t **errors_out,
                                size_t max_list);

/**
 * Free a quantum CSS code.
 */
void qcss_free(QuantumCSSCode *qcss);

/* ------------------------------------------------------------------ */
/*  L9.4 — Meta-Complexity: MCSP and List-Decoding                    */
/* ------------------------------------------------------------------ */

/**
 * The Minimum Circuit Size Problem (MCSP):
 * Given the truth table of a Boolean function f: {0,1}^m → {0,1}
 * of length N = 2^m, decide if f has a circuit of size ≤ s.
 *
 * Connection to list-decoding (Kabanets-Cai 2000):
 * If we can list-decode certain codes derived from circuit classes,
 * we can solve MCSP in subexponential time.
 *
 * This implementation encodes a circuit as a "codeword" and
 * uses list-decoding to find all small circuits consistent
 * with the truth table.
 */
typedef struct {
    size_t m;            /* number of input variables */
    size_t N;            /* truth table length = 2^m */
    size_t s;            /* circuit size bound */
    size_t num_gates;    /* number of gates in the candidate circuit */
    size_t *gate_types;  /* gate types (AND=0, OR=1, NOT=2) */
    size_t *gate_inputs; /* pairs of input wire indices */
    size_t *output_wire; /* which wire gives the output */
} CircuitCandidate;

/**
 * Encode a circuit as a truth table "codeword" by evaluating it
 * on all 2^m inputs.
 *
 * @param cc      Circuit description.
 * @param truth   Output truth table (length 2^m).
 * @return true if evaluation succeeded.
 */
bool mscp_encode_circuit(const CircuitCandidate *cc, bool *truth);

/**
 * MCSP via list-decoding: given a truth table T, find all circuits
 * of size ≤ s whose truth table agrees with T on at least (1-δ)N
 * positions.
 *
 * This is list-decoding the "circuit code" — the set of all
 * truth tables of circuits of size ≤ s.
 *
 * @param truth_in   Input truth table.
 * @param s          Circuit size bound.
 * @param delta      Error fraction tolerance.
 * @param circuits   Output candidate circuits.
 * @param max_list   Maximum number of candidates.
 * @return Number of candidates found.
 */
size_t mscp_list_decode_circuits(const bool *truth_in, size_t m, size_t s,
                                  double delta,
                                  CircuitCandidate **circuits,
                                  size_t max_list);

/**
 * Compute the MCSP list-decoding complexity estimate.
 *
 * Theorem (Informal): If AC0 circuits are list-decodable up to
 * error fraction 1/2-ε, then MCSP ∈ SUBEXP.
 *
 * @return Estimated time complexity exponent for MCSP.
 */
double mscp_complexity_estimate(size_t m, size_t s, double epsilon);

/* ------------------------------------------------------------------ */
/*  L9.5 — Locally Testable Codes (LTCs) and List-Decoding            */
/* ------------------------------------------------------------------ */

/**
 * A code is locally testable if there exists a randomized tester
 * that makes q = O(1) queries to a word and:
 *   - Accepts with prob 1 if word is a codeword
 *   - Rejects with prob Ω(δ(w,C)) if word is δ-far from the code
 *
 * The Hadamard code is both locally testable (3-query linearity test)
 * AND list-decodable (Goldreich-Levin).  This dual property is
 * central to many PCP constructions.
 */
typedef struct {
    size_t n;            /* block length */
    size_t q;            /* alphabet size */
    size_t queries;      /* number of queries for local test */
    double soundness;    /* rejection probability per unit distance */
    double list_radius;  /* list-decoding radius */
} LTCCode;

/**
 * Create an LTC descriptor for a Hadamard code.
 * Hadamard H_k: n=2^k, 3-query linearity test, soundness = 1/4.
 */
LTCCode *ltc_create_hadamard(size_t k);

/**
 * Local linearity test: query f(x), f(y), f(x⊕y) and check
 * f(x) ⊕ f(y) = f(x⊕y).  Reject if the check fails.
 *
 * @param oracle  Function to test.
 * @param k       Dimension (2^k = block length).
 * @param trials  Number of random trials.
 * @return Estimated distance to the Hadamard code.
 */
double ltc_test_linearity(bool (*oracle)(const bool *x), size_t k,
                           size_t trials);

/**
 * Combined LTC + list-decoding: first test locally, then
 * if the word passes the test, list-decode to find candidates.
 *
 * This is the paradigm behind the PCP theorem.
 *
 * @return Number of list-decoded candidates.
 */
size_t ltc_list_decode_pipeline(bool (*oracle)(const bool *x), size_t k,
                                 double epsilon, bool **decoded, size_t max_list);

/**
 * Free an LTC descriptor.
 */
void ltc_free(LTCCode *ltc);

#endif /* LD_RESEARCH_H */
