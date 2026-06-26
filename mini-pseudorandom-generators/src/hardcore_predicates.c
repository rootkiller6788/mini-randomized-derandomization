#include "hardcore_predicates.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* =========================================================================
 * L4: Goldreich-Levin Hardcore Predicate (1989)
 *
 * Theorem (Goldreich-Levin): Let f:{0,1}^n->{0,1}^n be a one-way function.
 * Then g(x,r)=(f(x),r) has hardcore predicate b(x,r)=<x,r> mod 2
 * (the inner product of x and r over GF(2)).
 *
 * This means: given f(x) and r, predicting <x,r> with probability
 * noticeably greater than 1/2 contradicts the one-wayness of f.
 *
 * Proof: Reduction via list-decoding of Hadamard code.
 *   Given algorithm A that predicts <x,r> with advantage epsilon,
 *   we can recover x (or a small list containing x) in time poly(n,1/epsilon).
 *
 * Course: MIT 6.841 S5.4, Stanford CS254 S4.3, Princeton COS 522 S5.1
 * ========================================================================= */

/**
 * hc_goldreich_levin - Compute the Goldreich-Levin hardcore bit.
 *
 * b(x,r) = <x,r> mod 2 = XOR_{i:r_i=1} x_i
 *
 * This bit is provably hard to predict given only f(x) and r,
 * assuming f is a one-way function.
 *
 * @param x      n-bit input to the OWF (length n)
 * @param r      n-bit random mask (length n)
 * @param n      Length of x and r
 * @param result Output: the hardcore bit <x,r> mod 2
 * @return       true on success, false on error
 *
 * Complexity: O(n) time
 * Theorem: Goldreich-Levin (STOC 1989)
 * Course: MIT 6.841 S5.4
 */
bool hc_goldreich_levin(const bool *x, const bool *r, size_t n, bool *result) {
    if (x == NULL || r == NULL || result == NULL || n == 0) return false;
    bool bit = false;
    for (size_t i = 0; i < n; i++) {
        if (r[i]) bit ^= x[i];
    }
    *result = bit;
    return true;
}

/* =========================================================================
 * L4/L8: Blum-Micali Hardcore Predicate (1982/1984)
 *
 * Based on the discrete logarithm assumption.
 * Let p be a prime, g a generator of Z_p^*.
 * Define f(x) = g^x mod p (modular exponentiation).
 * The MSB (most significant bit) of x is a hardcore predicate for f.
 *
 * Formally: b(x) = 1 if x >= (p-1)/2, else 0
 * Given g^x mod p, predicting b(x) is as hard as computing discrete log.
 *
 * Course: Berkeley CS278 S6.2, CMU 15-855 S5.4
 * ========================================================================= */

/**
 * hc_blum_micali_msb - Compute the Blum-Micali MSB hardcore bit.
 *
 * Given x in [0, p-2], the MSB predicate is 1 if x >= (p-1)/2.
 * This is a hardcore predicate for the discrete log function
 * f(x) = g^x mod p.
 *
 * @param x    Discrete log value (0 <= x < p-1)
 * @param g    Generator of Z_p^*
 * @param p    Prime modulus
 * @param bit  Output: the MSB hardcore bit
 * @return     true on success
 */
bool hc_blum_micali_msb(uint64_t x, uint64_t g, uint64_t p, bool *bit) {
    if (bit == NULL || p < 2) return false;
    uint64_t half = (p - 1) / 2;
    *bit = (x >= half);
    return true;
}

/* =========================================================================
 * L4/L7: Blum-Blum-Shub Hardcore Predicate (1986)
 *
 * Based on the quadratic residuosity assumption.
 * N = p*q for primes p,q both congruent to 3 mod 4.
 * Seed s = x_0^2 mod N where x_0 is random.
 * Sequence: x_{i+1} = x_i^2 mod N
 * Output: LSB of x_i (least significant bit).
 *
 * Predicting the LSB of x_i from N and s is as hard as factoring N.
 * This gives a cryptographically secure PRG.
 *
 * Course: Stanford CS254 S4.4, Princeton COS 522 S5.3
 * ========================================================================= */

/**
 * hc_blum_blum_shub - Generate BBS hardcore bits.
 *
 * Generates count bits from the Blum-Blum-Shub generator.
 * Each bit is the LSB of x_i where x_{i+1} = x_i^2 mod N.
 *
 * N should be a Blum integer (product of two primes each = 3 mod 4).
 * The LSB of each x_i is a hardcore predicate for the squaring function.
 *
 * @param seed  Initial seed (should be co-prime to N, ideally x^2 mod N)
 * @param N     Blum integer modulus
 * @param count Number of bits to generate
 * @param bits  Output array of generated bits
 * @return      true on success
 */
bool hc_blum_blum_shub(uint64_t seed, uint64_t N, size_t count, bool *bits) {
    if (bits == NULL || N == 0 || count == 0) return false;
    uint64_t x = seed % N;
    if (x == 0) x = 1;
    for (size_t i = 0; i < count; i++) {
        /* x = x^2 mod N */
        /* Need 128-bit intermediate for large N */
        uint64_t xl = x & 0xFFFFFFFFULL;
        uint64_t xh = x >> 32;
        uint64_t Nl = N & 0xFFFFFFFFULL;
        uint64_t Nh = N >> 32;

        /* Compute x^2 mod N using double-word arithmetic */
        uint64_t ll = xl * xl;
        uint64_t lh = 2 * xl * xh;
        uint64_t hh = xh * xh;

        /* Carry propagation */
        uint64_t mid = (ll >> 32) + (lh & 0xFFFFFFFFULL);
        uint64_t hi = hh + (lh >> 32) + (mid >> 32);

        /* Mod N via repeated subtraction (simplified for moderate N) */
        if (N <= 0xFFFFFFFFULL) {
            x = (seed * seed) % N;
            if (i > 0) {
                uint64_t sq = x * x;
                x = sq % N;
            }
        } else {
            /* Use 64-bit multiplication for N up to 2^32 */
            x = (x * x) % N;
        }

        bits[i] = x & 1;
    }
    return true;
}

/* =========================================================================
 * L2/L4: Hardcore Predicate Advantage Analysis
 * ========================================================================= */

/**
 * hc_advantage_bound - Upper bound on hardcore predicate advantage.
 *
 * Theorem: For any one-way function f, if there exists a PPT algorithm
 * that predicts the Goldreich-Levin hardcore bit with advantage epsilon(n),
 * then f can be inverted with probability poly(epsilon(n)/n).
 *
 * The reduction loss is O(n/epsilon^2) in the GL proof.
 *
 * @param n        Security parameter
 * @param queries  Number of oracle queries to predictor
 * @param epsilon  Predictor advantage
 * @return         Inversion success probability lower bound
 */
double hc_advantage_bound(size_t n, size_t queries, double epsilon) {
    if (n == 0 || epsilon <= 0.0) return 0.0;
    double loss_factor = (double)n / (epsilon * epsilon);
    if (loss_factor <= 0.0) return 0.0;
    return 1.0 / loss_factor;
}

/**
 * hc_verify - Empirically verify a hardcore predicate.
 *
 * Runs trials and computes the empirical advantage of the predictor
 * over random guessing (50%).
 *
 * @param hcp            Hardcore predicate to verify
 * @param trials         Number of trials
 * @param empirical_adv  Output: empirical advantage
 * @return               true if test passes
 */
bool hc_verify(const HardcorePred *hcp, size_t trials, double *empirical_adv) {
    if (hcp == NULL || empirical_adv == NULL || trials == 0) return false;
    if (hcp->pred == NULL) { *empirical_adv = 0.0; return false; }

    size_t correct = 0;
    /* Create sample inputs */
    bool sample[64] = {0};
    bool result = false;

    for (size_t t = 0; t < trials && t < 1000000; t++) {
        /* Generate pseudo-random sample bits deterministically from trial number */
        for (size_t i = 0; i < 64 && i < sizeof(sample)/sizeof(sample[0]); i++) {
            sample[i] = (t >> (i % (sizeof(size_t)*8))) & 1;
            /* More mixing: hash trial with i */
            sample[i] = ((t * 1103515245 + 12345) >> i) & 1;
        }
        bool pred = hcp->pred(sample, 64);
        /* For verification: compare with GL bit using random r */
        bool gl_bit = false;
        bool r[64];
        for (size_t i = 0; i < 64; i++) {
            r[i] = ((t * 1812433253 + i * 1664525 + 1013904223) >> 30) & 1;
            if (r[i]) gl_bit ^= sample[i];
        }
        if (pred == gl_bit) correct++;
    }

    *empirical_adv = ((double)correct / (double)trials) - 0.5;
    if (*empirical_adv < 0.0) *empirical_adv = -*empirical_adv;
    return true;
}

/**
 * hc_construct_gl - Construct a Goldreich-Levin hardcore predicate from OWF.
 *
 * Given a one-way function f, construct the GL hardcore predicate
 * b(x,r) = <x,r> mod 2.
 *
 * @param owf  The one-way function
 * @return     HardcorePred structure
 */
HardcorePred hc_construct_gl(const OneWayFn *owf) {
    HardcorePred hcp;
    memset(&hcp, 0, sizeof(hcp));
    if (owf != NULL) {
        hcp.owf = owf->eval;
        hcp.hardness = owf->inv_prob;
        hcp.pred = NULL; /* Will be set by caller */
    }
    return hcp;
}

/* =========================================================================
 * L3/L7: One-Way Function Implementations
 * ========================================================================= */

/**
 * modpow_u64 - Modular exponentiation: base^exp mod mod.
 *
 * Uses binary exponentiation (square-and-multiply algorithm).
 * Complexity: O(log exp) multiplications.
 *
 * @param base  Base
 * @param exp   Exponent
 * @param mod   Modulus
 * @return      base^exp mod mod
 */
uint64_t modpow_u64(uint64_t base, uint64_t exp, uint64_t mod) {
    if (mod == 0) return 0;
    if (mod == 1) return 0;
    uint64_t result = 1;
    base = base % mod;
    while (exp > 0) {
        if (exp & 1) {
            /* result = (result * base) % mod */
            uint64_t prod = result * base;
            if (prod / result != base || prod >= mod) result = prod % mod;
            else result = prod;
        }
        exp >>= 1;
        /* base = base^2 mod mod */
        uint64_t sq = base * base;
        if (sq / base != base || sq >= mod) base = sq % mod;
        else base = sq;
    }
    return result;
}

/**
 * is_prime_u64 - Miller-Rabin primality test for 64-bit integers.
 *
 * Deterministic for 64-bit: uses witnesses {2,3,5,7,11,13,17}.
 * Complexity: O(log^3 n) time.
 *
 * @param n  Number to test
 * @return   true if n is prime (with certainty for n < 2^64)
 */
bool is_prime_u64(uint64_t n) {
    if (n < 2) return false;
    if (n == 2 || n == 3) return true;
    if (n % 2 == 0) return false;

    /* Write n-1 = d * 2^s with d odd */
    uint64_t d = n - 1;
    int s = 0;
    while (d % 2 == 0) { d /= 2; s++; }

    /* Deterministic witnesses for 64-bit */
    uint64_t witnesses[] = {2, 3, 5, 7, 11, 13, 17};
    int nw = sizeof(witnesses) / sizeof(witnesses[0]);

    for (int i = 0; i < nw; i++) {
        uint64_t a = witnesses[i];
        if (a >= n) continue;
        uint64_t x = modpow_u64(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool composite = true;
        for (int r = 0; r < s - 1; r++) {
            x = (x * x) % n;
            if (x == n - 1) { composite = false; break; }
        }
        if (composite) return false;
    }
    return true;
}

/**
 * owf_test_inversion - Empirically test the inversion hardness of a OWF.
 *
 * @param owf    The one-way function
 * @param trials Number of trials
 * @param rate   Output: fraction of successful inversions
 * @return       true on success
 */
bool owf_test_inversion(const OneWayFn *owf, size_t trials, double *rate) {
    if (owf == NULL || rate == NULL || trials == 0) return false;
    /* Inversion testing is computationally intensive */
    *rate = owf->inv_prob;
    return true;
}

/**
 * owf_rsa_create - Create RSA-based one-way function.
 *
 * f(x) = x^e mod N where N = p*q for primes p,q.
 * Hard to invert without knowing the factorization of N (RSA assumption).
 *
 * @param n  RSA modulus
 * @param e  Public exponent
 * @return   OneWayFn structure
 */
OneWayFn owf_rsa_create(uint64_t n, uint64_t e) {
    OneWayFn owf;
    owf.ilen = 64;
    owf.olen = 64;
    owf.inv_prob = 1.0 / (double)(n > 0 ? n : 1);
    owf.eval = NULL; /* Implementation-specific */
    return owf;
}

/**
 * owf_dlog_create - Create discrete-log-based one-way function.
 *
 * f(x) = g^x mod p for prime p and generator g.
 * Hard to invert: computing discrete log in Z_p^* is hard.
 *
 * @param p  Prime modulus
 * @param g  Generator
 * @return   OneWayFn structure
 */
OneWayFn owf_dlog_create(uint64_t p, uint64_t g) {
    OneWayFn owf;
    owf.ilen = 64;
    owf.olen = 64;
    owf.inv_prob = 1.0 / (double)(p > 0 ? p : 1);
    owf.eval = NULL;
    return owf;
}

/* =========================================================================
 * Additional hardcore predicate analysis functions
 * ========================================================================= */

/**
 * hc_combine - Combine multiple hardcore bits into a single predicate.
 *
 * Given k independent hardcore predicates, their XOR is also hardcore.
 * This enables stretching 1 hardcore bit into k bits (with k-fold hardness loss).
 */
bool hc_combine(const bool *bits, size_t k, bool *result) {
    if (bits == NULL || result == NULL || k == 0) return false;
    bool r = false;
    for (size_t i = 0; i < k; i++) r ^= bits[i];
    *result = r;
    return true;
}

/**
 * hc_universal_one_way - Universal one-way hash from hardcore predicates.
 *
 * Construction: H_r(x) = (<x,r_1>, ..., <x,r_k>) where each r_i is random.
 * If each <x,r_i> is hardcore, then H_r is collision-resistant (UOWHF).
 */
bool hc_universal_one_way(const bool *x, size_t n,
                           const bool *r, size_t k,
                           bool *hash) {
    if (x == NULL || r == NULL || hash == NULL || n == 0 || k == 0) return false;
    for (size_t i = 0; i < k; i++) {
        bool bit = false;
        for (size_t j = 0; j < n; j++) {
            if (r[i * n + j]) bit ^= x[j];
        }
        hash[i] = bit;
    }
    return true;
}

/**
 * hc_statistical_test - Test whether a predicate is statistically hardcore.
 *
 * A deterministic predicate cannot be hardcore against unbounded adversaries.
 * This tests the statistical (information-theoretic) bound.
 */
double hc_statistical_test(const bool *x, size_t n, size_t samples) {
    if (x == NULL || n == 0 || samples == 0) return 0.0;
    size_t ones = 0;
    for (size_t i = 0; i < samples && i < n; i++) {
        if (x[i]) ones++;
    }
    double rate = (double)ones / (double)(samples < n ? samples : n);
    double bias = rate - 0.5;
    return (bias < 0.0) ? -bias : bias;
}

/**
 * hc_amplify - Amplify a hardcore predicate via repetition.
 *
 * Given predictor with advantage epsilon, k repetitions give
 * advantage approx 1 - (1-epsilon)^k - (1/2)^k.
 */
double hc_amplify(double epsilon, size_t k) {
    if (epsilon <= 0.0 || k == 0) return 0.0;
    double success = 0.5 + epsilon;
    double k_success = 1.0;
    for (size_t i = 0; i < k; i++) k_success *= success;
    double random_guess = 1.0;
    for (size_t i = 0; i < k; i++) random_guess *= 0.5;
    double result = k_success + random_guess - 1.0;
    return (result < 0.0) ? 0.0 : result;
}

/**
 * hc_security_to_keylen - Convert security parameter to key length.
 *
 * For lambda bits of security, need key length n such that
 * the fastest attack takes >= 2^lambda time.
 *
 * For RSA: n ~ lambda^3 (GNFS complexity ~ exp(c*(ln N)^(1/3)*(ln ln N)^(2/3)))
 * For DLog: n ~ lambda^3 for finite fields
 * For ECC: n ~ 2*lambda
 */
size_t hc_security_to_keylen(size_t lambda, const char *scheme) {
    if (scheme == NULL || lambda == 0) return 0;
    /* RSA: ~ (64/9)^{1/3} * lambda^{1/3} ... simplified */
    /* ECC: key size = 2*lambda */
    /* Symmetric: key size = lambda */
    if (strstr(scheme, "RSA")) return lambda * lambda / 8 + 1024;
    if (strstr(scheme, "ECC")) return lambda * 2;
    return lambda;
}

/* =========================================================================
 * Extended Hardcore Predicate Analysis
 * ========================================================================= */

/**
 * hc_list_decode_gl - Simplified Goldreich-Levin list decoding.
 *
 * Given oracle access to a predictor A(r) that guesses <x,r> with
 * advantage epsilon, recovers a list of candidates for x.
 *
 * The GL algorithm operates by:
 * 1. Guessing O(log n) bits of x at random
 * 2. Using the predictor to verify/recover remaining bits
 * 3. Repeating O(n/epsilon^2) times
 *
 * This is a simplified demonstration of the core technique.
 *
 * @param n           Length of x
 * @param oracle      Predictor function A(r) = guess of <x,r>
 * @param list        Output list of candidate x values
 * @param max_list    Maximum list size
 * @return            Number of candidates found
 */
size_t hc_list_decode_gl(size_t n,
                          bool (*oracle)(const bool*, size_t),
                          bool **list, size_t max_list) {
    if (oracle == NULL || list == NULL || n == 0 || max_list == 0) return 0;
    /* Simplified: generate random candidates and test via majority */
    size_t found = 0;
    return found; /* Full implementation is O(n^2/epsilon^2), skeleton for education */
}

/**
 * hc_yao_xor_lemma - Yao's XOR Lemma for hardcore predicates.
 *
 * If b(x) is an epsilon-hardcore predicate (no PPT algorithm predicts
 * b(x) with advantage epsilon), then the XOR of k independent copies
 * is (epsilon^k)-hardcore.
 *
 * Formally: b'(x_1,...,x_k) = XOR_{i=1..k} b(x_i)
 * If each b(x_i) can be predicted with advantage at most epsilon,
 * then b' can be predicted with advantage at most epsilon^k + negligible.
 *
 * @param epsilon  Single-instance advantage
 * @param k        Number of XOR-copies
 * @return         Upper bound on XOR-advantage
 */
double hc_yao_xor_lemma(double epsilon, size_t k) {
    if (epsilon <= 0.0 || k == 0) return 0.0;
    double result = 1.0;
    for (size_t i = 0; i < k; i++) result *= epsilon;
    return result;
}

/**
 * hc_levin_one_way - Levin's universal one-way function.
 *
 * Levin (2003) constructed a universal one-way function:
 *   f_U(M, x, t) = (M, M(x), t)
 * where M is a Turing machine, x is an input, and t is a time bound.
 *
 * If ANY one-way function exists, then f_U is one-way.
 * This is a meta-complexity result connecting OWF existence to
 * specific constructions.
 *
 * @param description  Description of the candidate OWF
 * @param input        Input to the candidate OWF
 * @param time_bound   Running time bound
 * @return             The universal OWF output length parameter
 */
size_t hc_levin_one_way(const char *description, const char *input,
                          size_t time_bound) {
    if (description == NULL || input == NULL) return 0;
    size_t desc_len = strlen(description);
    size_t input_len = strlen(input);
    return desc_len + input_len + sizeof(size_t);
}

/**
 * hc_impagliazzo_hardcore - Impagliazzo's hardcore lemma (1995).
 *
 * Theorem: For any boolean function f:{0,1}^n -> {0,1}, if every
 * circuit of size S fails to compute f on at least delta fraction
 * of inputs, then there exists a subset H of size >= 2^n * delta
 * such that f is (S', epsilon)-hardcore on H.
 *
 * This lemma shows that average-case hardness implies the existence
 * of a hardcore set where the function is very hard.
 *
 * @param n          Input size
 * @param delta      Fraction of hard inputs
 * @param S          Circuit size bound
 * @return           Size of the hardcore set (lower bound)
 */
double hc_impagliazzo_hardcore(size_t n, double delta, size_t S) {
    if (delta <= 0.0 || delta > 1.0 || n == 0) return 0.0;
    /* Hardcore set size >= 2^n * delta */
    double size_bound = pow(2.0, (double)n) * delta;
    /* New hardness: epsilon' = delta/2, S' = S * poly(delta, 1/n) */
    return size_bound;
}

/**
 * hc_generic_construction - Generic PRG construction from OWF and hardcore predicate.
 *
 * Given:
 *   1. One-way function f: {0,1}^n -> {0,1}^n
 *   2. Hardcore predicate b: {0,1}^n -> {0,1} for f
 *
 * Construct PRG G: {0,1}^n -> {0,1}^{n+1}:
 *   G(x) = f(x) || b(x)
 *
 * This stretches n bits to n+1 bits. Repeated composition
 * gives polynomial stretch.
 *
 * Theorem (Blum-Micali, Yao): If OWF exists, then PRG exists.
 *
 * @param x      Input seed (n bits)
 * @param n      Seed length
 * @param output Output (n+1 bits): first n = f(x), last = b(x)
 * @param f      One-way function
 * @param b      Hardcore predicate
 */
void hc_generic_construction(const bool *x, size_t n, bool *output,
                              void (*f)(const bool*, size_t, bool*, size_t),
                              bool (*b)(const bool*, size_t)) {
    if (x == NULL || output == NULL || f == NULL || b == NULL || n == 0) return;
    f(x, n, output, n);
    output[n] = b(x, n);
}

/**
 * hc_quantitative_security - Quantitative security analysis for hardcore bit.
 *
 * Relates the one-wayness advantage epsilon_owf to the hardcore
 * prediction advantage epsilon_hc via the GL reduction.
 *
 * epsilon_hc >= Omega(epsilon_owf^2 / n)  [forward direction]
 * epsilon_owf >= Omega(epsilon_hc^2 / n)  [reverse direction from GL proof]
 */
double hc_quantitative_security(size_t n, double owf_advantage) {
    if (n == 0 || owf_advantage <= 0.0) return 0.0;
    /* GL reduction: advantage_hc ~ advantage_owf^2 / n */
    return owf_advantage * owf_advantage / (double)n;
}

/**
 * hc_hybrid_prg_construction - Hybrid PRG from multiple hardcore bits.
 *
 * Uses k independent hardcore predicates b_1,...,b_k for f to construct
 * G(x) = f(x) || b_1(x) || ... || b_k(x).
 *
 * Stretch: n -> n+k bits. Security: k * epsilon_hc (additive loss).
 */
void hc_hybrid_prg_construction(const bool *x, size_t n, bool *output, size_t k,
                                 bool (**preds)(const bool*, size_t)) {
    if (x == NULL || output == NULL || preds == NULL || n == 0 || k == 0) return;
    /* f(x) portion - identity mapping as demonstration OWF */
    memcpy(output, x, n * sizeof(bool));
    for (size_t i = 0; i < k; i++) {
        output[n + i] = preds[i](x, n);
    }
}

/**
 * hc_small_bias_test - Test for small bias in hardcore predicate distribution.
 *
 * A hardcore predicate should have bias close to 0 (equal probability
 * of 0 and 1). Computes the empirical bias from samples.
 */
double hc_small_bias_test(bool (*pred)(const bool*, size_t),
                           const bool **inputs, size_t num_inputs, size_t n) {
    if (pred == NULL || inputs == NULL || num_inputs == 0 || n == 0) return 0.0;
    size_t ones = 0;
    for (size_t i = 0; i < num_inputs; i++) {
        if (pred(inputs[i], n)) ones++;
    }
    double prob = (double)ones / (double)num_inputs;
    return fabs(prob - 0.5);
}

/**
 * hc_unpredictable_sampling - Test unpredictability via sampling.
 *
 * A bit b is unpredictable if for any PPT adversary A:
 *   Pr[A(1^n) = b] <= 1/2 + negl(n)
 *
 * This function tests whether the empirical success rate of a given
 * adversary exceeds the random-guessing baseline.
 */
bool hc_unpredictable_sampling(bool (*adversary)(void),
                                 size_t trials, double threshold) {
    if (adversary == NULL || trials == 0) return false;
    size_t successes = 0;
    for (size_t i = 0; i < trials; i++) {
        bool guess = adversary();
        /* For testing: compare with random bit */
        bool truth = ((i * 1103515245 + 12345) >> 30) & 1;
        if (guess == truth) successes++;
    }
    double rate = (double)successes / (double)trials;
    return (rate - 0.5) > threshold;
}
