
/**
 * src/extractor_applications.c -- Applications of Randomness Extractors
 *
 * L7 Applications: BPP derandomization, cryptographic key derivation,
 *     privacy amplification, explicit extractor simulation
 * L8 Advanced: statistical test batteries, multi-source extraction
 * L9 Research Frontiers: quantum-proof extractors, network protocols
 *
 * Key applications:
 *   1. BPP derandomization via seed enumeration (Nisan-Wigderson, IW97)
 *   2. Cryptographic key derivation from weak Diffie-Hellman output
 *   3. Privacy amplification (Bennett-Brassard-Robert, 1988)
 *   4. Explicit extractor simulation with statistical validation
 *   5. Network extractor protocols (Kalai-Rao-Vadhan, 2009)
 *   6. Multi-source fault-tolerant extraction
 *   7. Quantum-proof extractor bounds
 *
 * Theorem citations:
 *   IW97: If E requires 2^{Omega(n)}-size circuits, then P = BPP
 *     (derandomization via NW PRG built from hard function)
 *   BBR88: Privacy amplification via universal hashing:
 *     Delta(Key | Eve's view, Uniform | Eve's view) <= 2^{-(k-m)/2}
 *   LHL: Cryptographic key derivation from any high-entropy source
 *   KRV09: Network extractors achieve multi-party randomness extraction
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>

#include "extractor_core.h"
#include "extractor_constructions.h"
#include "extractor_applications.h"

/*============================================================
 * CryptoKey helpers
 *============================================================*/

CryptoKey ck_create(size_t len) {
    CryptoKey ck;
    ck.len = len;
    ck.key = (uint8_t*)calloc(len, 1);
    return ck;
}

void ck_free(CryptoKey *ck) {
    if (ck && ck->key) {
        free(ck->key);
        ck->key = NULL;
        ck->len = 0;
    }
}

bool ck_eq(const CryptoKey *a, const CryptoKey *b) {
    if (!a || !b || a->len != b->len) return false;
    /* Constant-time comparison */
    uint8_t diff = 0;
    for (size_t i = 0; i < a->len; i++)
        diff |= a->key[i] ^ b->key[i];
    return diff == 0;
}

void ck_print(const CryptoKey *ck, const char *label) {
    if (!ck) return;
    if (label) printf("%s (%zu bytes): ", label, ck->len);
    for (size_t i = 0; i < ck->len && i < 64; i++)
        printf("%02x", ck->key[i]);
    if (ck->len > 64) printf("...");
    printf("\n");
}

/*============================================================
 * Section 1: BPP Derandomization
 *
 * Theorem (Nisan-Wigderson 1994, Impagliazzo-Wigderson 1997):
 *   If there exists a function f in E with circuit complexity
 *   2^{Omega(n)}, then P = BPP.
 *
 * Derandomization procedure:
 *   1. Build PRG: {0,1}^{O(log n)} -> {0,1}^{poly(n)} from f
 *   2. For each deterministic input x, enumerate all PRG seeds
 *   3. Run BPP algorithm on each PRG output
 *   4. Take majority vote
 *
 * Complexity: 2^{O(log n)} = poly(n) time.
 *
 * This implementation: enumeration-based (feasible for small n).
 *============================================================*/

bool ext_derandomize_bpp(bool (*alg)(const char*, size_t,
                                      const bool*, size_t),
                          const char *inp, size_t ilen,
                          size_t bits, bool *result) {
    if (!alg || !inp || !result || bits == 0) return false;
    if (bits > 24) {
        fprintf(stderr, "ext_derandomize_bpp: %zu bits too many (>24)\n", bits);
        return false;
    }
    size_t nseeds = (size_t)1 << bits;
    size_t ones = 0, zeros = 0;
    bool *seed = (bool*)malloc(bits * sizeof(bool));
    if (!seed) return false;

    for (size_t s = 0; s < nseeds; s++) {
        /* Generate seed from counter */
        for (size_t i = 0; i < bits; i++)
            seed[i] = (s >> i) & 1;
        if (alg(inp, ilen, seed, bits))
            ones++;
        else
            zeros++;
    }
    free(seed);

    *result = (ones > zeros);
    size_t decisive = (ones > zeros) ? ones : zeros;
    return ((double)decisive / (double)nseeds >= 0.666);
}

/**
 * NW-style derandomization: use a shorter seed via PRG.
 *
 * Instead of enumerating all 2^n random bits, we:
 *   1. Generate a PRG: {0,1}^s -> {0,1}^n where s << n
 *   2. Enumerate only 2^s seeds
 *
 * For this demo, the PRG is the NW generator from extractor_constructions.
 */
bool ext_bpp_nw_derandomize(bool (*alg)(const char*, size_t,
                                         const bool*, size_t),
                             const char *inp, size_t ilen,
                             size_t seed_len, size_t prg_out_len,
                             bool *result, uint64_t *rng_seed) {
    if (!alg || !inp || !result || seed_len == 0 || prg_out_len == 0)
        return false;
    if (seed_len > 20) {
        fprintf(stderr, "ext_bpp_nw_derandomize: seed_len=%zu too large\n",
                seed_len);
        return false;
    }
    size_t nseeds = (size_t)1 << seed_len;
    size_t ones = 0, zeros = 0;
    bool *seed = (bool*)malloc(seed_len * sizeof(bool));
    bool *prg_out = (bool*)malloc(prg_out_len * sizeof(bool));
    if (!seed || !prg_out) { free(seed); free(prg_out); return false; }

    for (size_t s = 0; s < nseeds; s++) {
        for (size_t i = 0; i < seed_len; i++)
            seed[i] = (s >> i) & 1;

        /* Use LCG to generate prg_out from seed (simple PRG for demo) */
        uint64_t prg_state = *rng_seed ^ s;
        for (size_t i = 0; i < prg_out_len; i++) {
            prg_state = prg_state * 6364136223846793005ULL
                        + 1442695040888963407ULL;
            prg_out[i] = (prg_state >> 63) & 1;
        }

        if (alg(inp, ilen, prg_out, prg_out_len))
            ones++;
        else
            zeros++;
    }
    free(seed); free(prg_out);

    *result = (ones > zeros);
    size_t decisive = (ones > zeros) ? ones : zeros;
    return ((double)decisive / (double)nseeds >= 0.666);
}

/*============================================================
 * Section 2: Cryptographic Key Derivation
 *
 * Problem: Alice and Bob perform Diffie-Hellman and get
 *   g^{ab} mod p, which has high min-entropy but is not
 *   uniformly random. They want a uniform cryptographic key.
 *
 * Solution: Use a randomness extractor.
 *   key = Ext(g^{ab} mod p, seed)
 * where seed can be public (sent in the clear).
 *
 * Theorem (LHL, ILL 1989):
 *   If H_inf(weak_secret) >= k bits, then the extracted key
 *   is 2^{-(k - key_len*8)/2}-close to uniform.
 *
 * Implementation: use Carter-Wegman universal hash as extractor.
 *============================================================*/

CryptoKey ext_key_derivation(const uint8_t *weak_secret, size_t slen,
                              size_t key_len) {
    CryptoKey ck = ck_create(key_len);
    if (!weak_secret || slen == 0) return ck;

    /* Interpret weak_secret bytes as a large integer for hashing */
    size_t n_bits = slen * 8;
    size_t m_bits = key_len * 8;
    if (n_bits > 64) n_bits = 64; /* Limit for this demo */
    if (m_bits > 64) m_bits = 64;

    UniversalHash uh = uh_create_carter_wegman(n_bits, m_bits);

    /* Convert bytes to integer */
    size_t x = 0;
    for (size_t i = 0; i < slen && i < 8; i++)
        x |= ((size_t)weak_secret[i]) << (i * 8);

    /* Use a fixed seed for determinism (in practice, use random seed) */
    size_t a_idx = 42 % uh.num_a;
    size_t b_idx = 137 % uh.num_b;
    size_t hash_val = uh_eval(&uh, a_idx, b_idx, x);

    /* Write hash output to key bytes */
    for (size_t i = 0; i < key_len; i++)
        ck.key[i] = (uint8_t)((hash_val >> (i * 8)) & 0xFF);

    uh_free(&uh);
    return ck;
}

double ext_key_derivation_bound(size_t slen_bits, double min_entropy_bits,
                                 size_t key_len_bits) {
    (void)slen_bits;
    if (min_entropy_bits <= (double)key_len_bits) return 1.0;
    return pow(2.0, -(min_entropy_bits - (double)key_len_bits) / 2.0);
}

/*============================================================
 * Section 3: Privacy Amplification (Bennett-Brassard-Robert 1988)
 *
 * Scenario (QKD / information-theoretic key agreement):
 *   - Alice and Bob share a string w (partially known to Eve)
 *   - H_inf(w | Eve's information) >= k
 *   - They publicly choose a random seed s
 *   - Both compute key = Ext(w, s)
 *   - Result: Delta((s, key), (s, uniform)) <= 2^{-(k-m)/2}
 *
 * This means Eve's information about the key is exponentially
 * small in the entropy gap.
 *============================================================*/

bool ext_privacy_amplification(const uint8_t *shared, size_t slen,
                                const uint8_t *seed, size_t seed_len,
                                uint8_t *key, size_t key_len) {
    if (!shared || !seed || !key || slen == 0 || key_len == 0)
        return false;

    /* Combine shared secret and seed via XOR-then-extract */
    size_t n_bits = slen * 8;
    size_t m_bits = key_len * 8;
    if (n_bits > 64) n_bits = 64;
    if (m_bits > 64) m_bits = 64;

    UniversalHash uh = uh_create_carter_wegman(n_bits, m_bits);

    /* Convert shared bytes to integer */
    size_t x = 0;
    for (size_t i = 0; i < slen && i < 8; i++)
        x |= ((size_t)shared[i]) << (i * 8);

    /* Derive hash parameters from seed */
    size_t seed_val = 0;
    for (size_t i = 0; i < seed_len && i < 8; i++)
        seed_val |= ((size_t)seed[i]) << (i * 8);

    size_t a_idx = seed_val % uh.num_a;
    size_t b_idx = (seed_val / uh.num_a) % uh.num_b;
    size_t hash_val = uh_eval(&uh, a_idx, b_idx, x);

    /* Write key bytes */
    for (size_t i = 0; i < key_len; i++)
        key[i] = (uint8_t)((hash_val >> (i * 8)) & 0xFF);

    uh_free(&uh);
    return true;
}

double ext_privacy_amplification_bound(double k_bits, size_t m_bits) {
    if (k_bits <= (double)m_bits) return 1.0;
    return pow(2.0, -(k_bits - (double)m_bits) / 2.0);
}

/*============================================================
 * Section 4: Explicit Extractor Simulation
 *
 * Empirically test an extractor's output uniformity by:
 *   1. Generate random k-sources
 *   2. Run extractor
 *   3. Apply statistical tests (chi-squared, frequency, runs)
 *   4. Compare observed statistical distance to theoretical bound
 *============================================================*/

double ext_simulate_explicit_extractor(size_t n, size_t k, size_t d,
                                        size_t m, double eps,
                                        size_t trials) {
    if (n == 0 || n > 12 || k > n || m == 0 || trials == 0) return -1.0;
    (void)d; (void)eps;

    size_t nvals = (size_t)1 << m;
    if (nvals > 65536) nvals = 65536;
    size_t *freq = (size_t*)calloc(nvals, sizeof(size_t));
    if (!freq) return -1.0;

    bool *src = (bool*)malloc(n * sizeof(bool));
    bool *seed = (bool*)malloc(d * sizeof(bool));
    bool *out = (bool*)malloc(m * sizeof(bool));
    if (!src || !seed || !out) {
        free(freq); free(src); free(seed); free(out); return -1.0;
    }

    uint64_t rng = 42;
    Extractor e = ext_create(n, d, m, eps);
    if (!ext_is_valid(&e)) {
        free(freq); free(src); free(seed); free(out); return -1.0;
    }

    for (size_t t = 0; t < trials; t++) {
        /* Generate source with min-entropy approx k */
        double bias = 1.0 - pow(2.0, -(double)k / (double)n);
        for (size_t i = 0; i < n; i++)
            src[i] = ((double)(ext_lcg_next(&rng) >> 11) / (1ULL << 53)) < bias;
        /* Uniform seed */
        for (size_t i = 0; i < d; i++)
            seed[i] = (ext_lcg_next(&rng) >> 63) & 1;
        if (ext_eval(&e, src, seed, out)) {
            size_t val = 0;
            for (size_t i = 0; i < m; i++)
                if (out[i]) val |= ((size_t)1) << i;
            if (val < nvals) freq[val]++;
        }
    }
    free(src); free(seed); free(out);

    /* Estimate statistical distance to uniform */
    size_t total = 0;
    for (size_t i = 0; i < nvals; i++) total += freq[i];
    if (total == 0) { free(freq); return -1.0; }

    double delta = 0.0;
    double uniform_prob = 1.0 / (double)nvals;
    for (size_t i = 0; i < nvals; i++) {
        double observed = (double)freq[i] / (double)total;
        delta += fabs(observed - uniform_prob);
    }
    delta *= 0.5;
    free(freq);
    return delta;
}

double ext_statistical_test_battery(const bool *output, size_t n_samples,
                                     size_t m) {
    if (!output || n_samples == 0 || m == 0) return -1.0;
    (void)output; (void)n_samples; (void)m;

    /* Monobit (frequency) test: count ones */
    size_t ones = 0;
    for (size_t i = 0; i < n_samples * m; i++)
        if (output[i]) ones++;
    double p_monobit = erfc(fabs((double)ones - (double)(n_samples*m)/2.0)
                            / sqrt((double)(n_samples*m) / 2.0));

    /* Runs test: count runs of consecutive identical bits */
    size_t runs = 1;
    for (size_t i = 1; i < n_samples * m; i++)
        if (output[i] != output[i-1]) runs++;
    double pi = (double)ones / (double)(n_samples * m);
    double runs_expected = 1.0 + 2.0 * (double)(n_samples*m-1) * pi * (1.0 - pi);
    double runs_std = sqrt(2.0 * (double)(n_samples*m) * pi * (1.0-pi)
                           * (1.0 - 3.0*pi*(1.0-pi)));
    double p_runs = erfc(fabs((double)runs - runs_expected) / runs_std);

    /* Chi-squared test on m-bit blocks */
    size_t nvals = (size_t)1 << m;
    if (nvals > n_samples) nvals = n_samples;
    size_t *freq = (size_t*)calloc(nvals, sizeof(size_t));
    if (freq) {
        for (size_t s = 0; s < n_samples; s++) {
            size_t val = 0;
            for (size_t i = 0; i < m; i++)
                if (output[s * m + i]) val |= ((size_t)1) << i;
            if (val < nvals) freq[val]++;
        }
        double expected = (double)n_samples / (double)nvals;
        double chi2 = 0.0;
        for (size_t i = 0; i < nvals; i++) {
            double diff = (double)freq[i] - expected;
            chi2 += diff * diff / expected;
        }
        /* Approximate p-value from chi-squared distribution */
        /* Using Wilson-Hilferty transformation */
        double z = (pow(chi2 / (double)(nvals - 1), 1.0/3.0)
                    - (1.0 - 2.0/(9.0*(double)(nvals-1))))
                   / sqrt(2.0 / (9.0 * (double)(nvals - 1)));
        double p_chi2 = erfc(fabs(z));
        free(freq);
        /* Bonferroni: return min p-value * 3 */
        double p_min = p_monobit;
        if (p_runs < p_min) p_min = p_runs;
        if (p_chi2 < p_min) p_min = p_chi2;
        return p_min * 3.0;
    }

    /* Return min of available p-values */
    double p_min = p_monobit;
    if (p_runs < p_min) p_min = p_runs;
    return p_min;
}

/*============================================================
 * Section 5: Network Extractor Protocol
 *
 * Model (KRV 2009):
 *   - t parties, each with a local weak source X_i
 *   - Adversary controls up to f parties
 *   - Honest parties want to extract a common random string
 *
 * Protocol (star topology):
 *   1. Each party sends local source to central node
 *   2. Central node runs extractor on concatenated sources
 *   3. Broadcasts result
 *
 * Fault tolerance: even if f sources are bad, as long as
 * t-f sources have min-entropy >= k, extraction succeeds.
 *============================================================*/

bool ext_network_extractor_simulate(size_t n, size_t t, size_t d,
                                     size_t m, size_t trials) {
    if (n == 0 || t == 0 || m == 0 || trials == 0) return false;
    (void)d;

    bool *concatenated = (bool*)malloc(n * t * sizeof(bool));
    bool *seed = (bool*)malloc(d * sizeof(bool));
    bool *out = (bool*)malloc(m * sizeof(bool));
    if (!concatenated || !seed || !out) {
        free(concatenated); free(seed); free(out); return false;
    }

    uint64_t rng = 99;
    size_t successes = 0;
    Extractor e = ext_create(n * t, d, m, 0.1);

    for (size_t tr = 0; tr < trials; tr++) {
        /* Generate t independent sources */
        for (size_t i = 0; i < t; i++) {
            double bias = 0.3 + 0.4 * ((double)(ext_lcg_next(&rng) & 0xFFFF) / 0x10000);
            for (size_t j = 0; j < n; j++)
                concatenated[i * n + j] =
                    ((double)(ext_lcg_next(&rng)>>11)/(1ULL<<53)) < bias;
        }
        for (size_t i = 0; i < d; i++)
            seed[i] = (ext_lcg_next(&rng) >> 63) & 1;

        if (ext_eval(&e, concatenated, seed, out))
            successes++;
    }

    free(concatenated); free(seed); free(out);
    return ((double)successes / (double)trials) >= 0.95;
}

bool ext_multi_source_fault_tolerant(const bool **sources, size_t n,
                                      size_t t, size_t f, size_t d,
                                      size_t m, bool *output) {
    if (!sources || !output || n == 0 || t == 0 || m == 0) return false;
    if (f >= t) return false; /* All parties corrupted */

    /* Concatenate all sources */
    bool *concat = (bool*)malloc(n * t * sizeof(bool));
    bool *seed = (bool*)malloc(d * sizeof(bool));
    if (!concat || !seed) { free(concat); free(seed); return false; }

    for (size_t i = 0; i < t; i++)
        memcpy(concat + i * n, sources[i], n * sizeof(bool));

    /* Generate public seed (fixed for demo) */
    for (size_t i = 0; i < d; i++) seed[i] = 0;

    Extractor e = ext_create(n * t, d, m, 0.1);
    bool ok = ext_eval(&e, concat, seed, output);

    free(concat); free(seed);
    return ok;
}

double ext_quantum_proof_bound(double k_bits, size_t m_bits) {
    /* Quantum-proof extraction (Tomamichel et al. 2011):
     * When the adversary holds quantum side information E,
     * the extractor with universal hashing still works.
     *
     * Delta( rho_{KEY,E}, rho_{U,E} ) <= 2 * 2^{-(k-m)/4}
     * (slightly worse than classical but qualitatively same).
     */
    if (k_bits <= (double)m_bits) return 1.0;
    return 2.0 * pow(2.0, -(k_bits - (double)m_bits) / 4.0);
}

/*============================================================
 * Section 6: Self-test
 *============================================================*/

int ext_apps_self_test(void) {
    int f = 0;

    /* CryptoKey create/free */
    {
        CryptoKey ck = ck_create(16);
        if (ck.len != 16 || !ck.key) f++;
        for (size_t i = 0; i < 16; i++)
            if (ck.key[i] != 0) f++;
        ck_free(&ck);
    }

    /* CryptoKey equality */
    {
        CryptoKey a = ck_create(4), b = ck_create(4);
        a.key[0] = 0xAB; b.key[0] = 0xAB;
        a.key[1] = 0xCD; b.key[1] = 0xCD;
        if (!ck_eq(&a, &b)) f++;
        b.key[1] = 0xCE;
        if (ck_eq(&a, &b)) f++;
        ck_free(&a); ck_free(&b);
    }

    /* Key derivation */
    {
        uint8_t weak[8] = {0x12,0x34,0x56,0x78,0x9A,0xBC,0xDE,0xF0};
        CryptoKey key = ext_key_derivation(weak, 8, 4);
        if (key.len != 4) f++;
        ck_free(&key);
    }

    /* Privacy amplification */
    {
        uint8_t shared[8] = {0xFF,0xEE,0xDD,0xCC,0xBB,0xAA,0x99,0x88};
        uint8_t seed[4] = {0xCA,0xFE,0xBA,0xBE};
        uint8_t key[4] = {0};
        if (!ext_privacy_amplification(shared, 8, seed, 4, key, 4)) f++;
    }

    /* Privacy amplification bound */
    {
        double eps = ext_privacy_amplification_bound(128.0, 64);
        /* eps = 2^{-32} ~ 2.3e-10 */
        if (eps < 1e-15 || eps > 1e-5) f++;
    }

    /* Key derivation bound */
    {
        double eps = ext_key_derivation_bound(256, 128.0, 64);
        if (eps < 1e-15 || eps > 1e-5) f++;
    }

    /* Quantum proof bound */
    {
        double eps = ext_quantum_proof_bound(128.0, 64);
        if (eps < 1e-12 || eps > 1e-1) f++;
    }

    return f;
}
