#ifndef EXT_APPS_H
#define EXT_APPS_H

/**
 * mini-randomness-extractors — Applications
 *
 * L7 Applications: BPP derandomization, cryptographic key derivation,
 *                   privacy amplification, explicit extractor simulation
 * L9 Research Frontiers: network extractors, multi-source extraction
 *
 * References:
 *   Nisan & Wigderson (1994) "Hardness vs Randomness"
 *   Bennett, Brassard, Robert (1988) "Privacy Amplification"
 *   Impagliazzo & Wigderson (1997) "P = BPP if E requires
 *     exponential circuits"
 *   Kalai, Rao, Vadhan (2009) "Network Extractor Protocols"
 */

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>

/*──────────────────────────────────────────────────────────────────────
 * L7: Application types
 *──────────────────────────────────────────────────────────────────────*/

/** Entropy source specification for applications.
 */
typedef struct {
    size_t n;        /**< Source length in bits */
    double k;        /**< Guaranteed min-entropy in bits */
} EntropySource;

/** Cryptographic key material.
 */
typedef struct {
    uint8_t *key;    /**< Key bytes (heap-allocated) */
    size_t len;      /**< Key length in bytes */
} CryptoKey;

/** Allocate and return a CryptoKey of given length (filled with zeros). */
CryptoKey ck_create(size_t len);
/** Free a CryptoKey. */
void ck_free(CryptoKey *ck);
/** Compare two keys in constant time (for security). Returns true if equal. */
bool ck_eq(const CryptoKey *a, const CryptoKey *b);
/** Hex dump a key to stdout. */
void ck_print(const CryptoKey *ck, const char *label);

/*──────────────────────────────────────────────────────────────────────
 * L7: BPP derandomization
 *──────────────────────────────────────────────────────────────────────*/

/** Derandomize a BPP algorithm by enumerating all seeds.
 *
 *  Theorem (NW, IW): If there exists a function in E with circuit
 *  complexity 2^{Ω(n)}, then P = BPP.
 *
 *  This function simulates the derandomization: given a BPP algorithm
 *  A(x, r) that uses n random bits and errs with prob < 1/3, we run
 *  A on all possible seeds r (enumerating the output of the extractor)
 *  and take majority vote. For small n this is feasible.
 *
 *  alg:    the BPP algorithm A(input, input_len, random_seed, seed_len) → bool
 *  inp:    the deterministic input
 *  ilen:   length of input
 *  bits:   number of random bits A expects
 *  result: output majority vote
 *
 *  Returns true if majority is decisive (≥ 2/3 agreement).
 */
bool ext_derandomize_bpp(bool (*alg)(const char*, size_t,
                                      const bool*, size_t),
                          const char *inp, size_t ilen,
                          size_t bits, bool *result);

/** Nisan-Wigderson "Hardness vs Randomness" derandomization.
 *  Uses a PRG to stretch a short truly random seed into many
 *  pseudorandom bits, which are fed to the BPP algorithm.
 *
 *  seed_len: length of truly random seed used by PRG
 *  prg_out_len: number of pseudorandom bits generated
 */
bool ext_bpp_nw_derandomize(bool (*alg)(const char*, size_t,
                                         const bool*, size_t),
                             const char *inp, size_t ilen,
                             size_t seed_len, size_t prg_out_len,
                             bool *result, uint64_t *rng_seed);

/*──────────────────────────────────────────────────────────────────────
 * L7: Cryptographic key derivation
 *──────────────────────────────────────────────────────────────────────*/

/** Derive a cryptographic key from a weak secret using a
 *  randomness extractor (Leftover Hash Lemma construction).
 *
 *  weak_secret: input with min-entropy ≥ k (e.g., Diffie-Hellman output)
 *  slen: length of weak_secret in bytes
 *  key_len: desired output key length in bytes
 *
 *  Uses Carter-Wegman universal hash as the extractor.
 *  The seed is drawn from system entropy; for deterministic testing,
 *  a fixed seed is used.
 *
 *  Theorem (LHL): If H∞(weak_secret) ≥ k, then the derived key is
 *  ≤ 2^{-(k_bits - 8*key_len)/2}-close to uniform.
 */
CryptoKey ext_key_derivation(const uint8_t *weak_secret, size_t slen,
                              size_t key_len);

/** Computes the leftover hash bound for given parameters.
 *  Returns the statistical distance guarantee (ε).
 */
double ext_key_derivation_bound(size_t slen_bits, double min_entropy_bits,
                                 size_t key_len_bits);

/*──────────────────────────────────────────────────────────────────────
 * L7: Privacy amplification (Bennett-Brassard-Robert)
 *──────────────────────────────────────────────────────────────────────*/

/** Privacy amplification: Alice and Bob share a partially secret
 *  string w (Eve has side information limiting H∞(w|Eve) ≥ k).
 *  They use a public random seed s to extract a shorter,
 *  almost-uniform key that Eve knows almost nothing about.
 *
 *  Protocol:
 *    1. Alice and Bob share weak_secret w (slen bytes)
 *    2. They publicly agree on seed (seed_len bytes)
 *    3. Both compute key = Ext(w, seed) via universal hash
 *    4. Security: Δ((Seed, Key), (Seed, Uniform)) ≤ 2^{-(k_bits - 8*key_len)/2}
 *
 *  shared: the weak shared secret
 *  slen: length of shared in bytes
 *  seed: public random seed
 *  seed_len: seed length in bytes
 *  key: output buffer (must be key_len bytes)
 *  key_len: desired derived key length
 *
 *  Returns true on success, false if parameters invalid.
 */
bool ext_privacy_amplification(const uint8_t *shared, size_t slen,
                                const uint8_t *seed, size_t seed_len,
                                uint8_t *key, size_t key_len);

/** Compute the privacy amplification security bound.
 *  Given shared secret min-entropy k, seed length d, output length m:
 *  ε ≤ 2^{-(k - m)/2}  (from LHL).
 */
double ext_privacy_amplification_bound(double k_bits, size_t m_bits);

/*──────────────────────────────────────────────────────────────────────
 * L7: Explicit extractor simulation
 *──────────────────────────────────────────────────────────────────────*/

/** Simulate an explicit extractor and empirically measure its
 *  output uniformity.
 *
 *  Generates random k-sources, runs the extractor, and estimates
 *  the statistical distance to uniform via chi-squared testing.
 *
 *  n: source length (must be ≤ 12 for exhaustive enumeration)
 *  k: min-entropy
 *  d: seed length
 *  m: output length
 *  eps: theoretical error bound
 *  trials: number of random trials
 *
 *  Returns estimated statistical distance.
 */
double ext_simulate_explicit_extractor(size_t n, size_t k, size_t d,
                                        size_t m, double eps,
                                        size_t trials);

/** Run a battery of statistical tests on extractor output.
 *  Tests: frequency (monobit), runs, chi-squared, serial correlation.
 *  Returns the minimum p-value across all tests (Bonferroni-adjusted).
 */
double ext_statistical_test_battery(const bool *output, size_t n_samples,
                                     size_t m);

/*──────────────────────────────────────────────────────────────────────
 * L9: Network extractor protocols
 *──────────────────────────────────────────────────────────────────────*/

/** Simulate a network extractor protocol (Kalai-Rao-Vadhan).
 *
 *  In a network of t parties, each with a local weak source,
 *  they run a protocol to extract a common random string.
 *
 *  Model:
 *    - t parties, each i has source X_i with H∞(X_i) ≥ k
 *    - Allowed communication rounds within a network topology
 *    - Goal: extract m uniform bits shared by all honest parties
 *
 *  This simulation uses a star-topology protocol: all parties send
 *  their source to a central party who runs an extractor.
 *
 *  n: source length per party
 *  t: number of parties
 *  d: seed length (public)
 *  m: output length per party
 *  trials: Monte Carlo trials
 *
 *  Returns true if the protocol succeeds within ε tolerance.
 */
bool ext_network_extractor_simulate(size_t n, size_t t, size_t d,
                                     size_t m, size_t trials);

/** Multi-source extractor (fault-tolerant).
 *  Even if up to f parties are corrupted, extract from the
 *  remaining good sources.
 */
bool ext_multi_source_fault_tolerant(const bool **sources, size_t n,
                                      size_t t, size_t f, size_t d,
                                      size_t m, bool *output);

/** Simulate a quantum-proof extractor.
 *  In the presence of quantum side information, the extractor must
 *  produce bits that look uniform even to a quantum adversary.
 *  This is guaranteed by the same LHL when using universal hashing.
 *  This function computes the quantum-proof security bound.
 */
double ext_quantum_proof_bound(double k_bits, size_t m_bits);

#endif /* EXT_APPS_H */
