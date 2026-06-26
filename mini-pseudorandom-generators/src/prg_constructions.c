#include "prg_constructions.h"
#include "next_bit_tests.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* =========================================================================
 * L1/L5: Linear Congruential Generator (LCG)
 * X_{n+1} = (a * X_n + c) mod m. Hull-Dobell Theorem for full period.
 * ========================================================================= */

LCG lcg_create(uint64_t seed, uint64_t a, uint64_t c, uint64_t m) {
    LCG lcg; lcg.a = a; lcg.c = c; lcg.m = m;
    lcg.state = (m > 0) ? seed % m : seed;
    if (lcg.state == 0 && c == 0) lcg.state = 1;
    return lcg;
}

uint64_t lcg_next(LCG *lcg) {
    if (lcg == NULL) return 0;
    uint64_t a = lcg->a, s = lcg->state, c = lcg->c, m = lcg->m;
    /* Compute (a*s + c) with overflow handling */
    if (m > 0 && m <= 0xFFFFFFFFULL) {
        lcg->state = ((a % m) * (s % m) + c) % m;
    } else if (m > 0) {
        lcg->state = a * s + c;
        if (lcg->state >= m) lcg->state %= m;
    } else {
        lcg->state = a * s + c;
    }
    return lcg->state;
}

double lcg_quality_test(const LCG *lcg, size_t n) {
    if (lcg == NULL || n == 0) return 0.0;
    LCG copy = *lcg;
    double sum = 0.0, sum_sq = 0.0;
    for (size_t i = 0; i < n; i++) {
        double val = (double)lcg_next(&copy);
        sum += val; sum_sq += val * val;
    }
    double mean = sum / (double)n;
    double expected = (lcg->m == 0) ? 9.22e18 : (double)lcg->m / 2.0;
    double dev = expected > 0 ? fabs(mean - expected) / expected : 1.0;
    double score = 1.0 - dev;
    return (score < 0.0) ? 0.0 : score;
}

/* =========================================================================
 * L5: Mersenne Twister MT19937 (Matsumoto-Nishimura 1998)
 * Period 2^{19937}-1. State: 624x32-bit. Twisted GFSR over GF(2).
 * ========================================================================= */

#define MT_N 624
#define MT_M 397
#define MT_MATRIX_A 0x9908B0DFUL
#define MT_UPPER 0x80000000UL
#define MT_LOWER 0x7FFFFFFFUL

MT19937 mt19937_create(uint64_t seed) {
    MT19937 mt; mt.seed = seed; mt.index = MT_N;
    mt.state[0] = seed & 0xFFFFFFFFUL;
    for (int i = 1; i < MT_N; i++) {
        mt.state[i] = (1812433253UL * (mt.state[i-1] ^ (mt.state[i-1] >> 30)) + i);
        mt.state[i] &= 0xFFFFFFFFUL;
    }
    return mt;
}

static void mt_twist(MT19937 *mt) {
    static const uint64_t mag[2] = {0UL, MT_MATRIX_A};
    for (int i = 0; i < MT_N; i++) {
        uint64_t x = (mt->state[i] & MT_UPPER) | (mt->state[(i+1)%MT_N] & MT_LOWER);
        mt->state[i] = mt->state[(i+MT_M)%MT_N] ^ (x>>1) ^ mag[x&1UL];
    }
    mt->index = 0;
}

uint64_t mt19937_extract(MT19937 *mt) {
    if (mt == NULL) return 0;
    if (mt->index >= MT_N) mt_twist(mt);
    uint64_t y = mt->state[mt->index++];
    y ^= (y>>11); y ^= (y<<7) & 0x9D2C5680UL;
    y ^= (y<<15) & 0xEFC60000UL; y ^= (y>>18);
    return y & 0xFFFFFFFFUL;
}

void mt19937_generate_bits(MT19937 *mt, bool *output, size_t n) {
    if (mt == NULL || output == NULL || n == 0) return;
    for (size_t i = 0; i < n; i++) {
        if (i % 32 == 0) {
            uint64_t val = mt19937_extract(mt);
            for (size_t j = 0; j < 32 && (i+j) < n; j++)
                output[i+j] = (val >> j) & 1;
        }
    }
}

/* =========================================================================
 * L3/L5: Xorshift128 (Marsaglia 2003). Period 2^{128}-1.
 * ========================================================================= */

Xorshift128 xorshift128_create(uint64_t seed1, uint64_t seed2) {
    Xorshift128 xs;
    xs.s[0] = (seed1==0) ? 123456789ULL : seed1;
    xs.s[1] = (seed2==0) ? 987654321ULL : seed2;
    return xs;
}

uint64_t xorshift128_next(Xorshift128 *xs) {
    if (xs == NULL) return 0;
    uint64_t s1 = xs->s[0], s0 = xs->s[1];
    uint64_t result = s0 + s1;
    xs->s[0] = s0;
    s1 ^= s1 << 23; s1 ^= s1 >> 17;
    s1 ^= s0; s1 ^= s0 >> 26;
    xs->s[1] = s1;
    return result;
}

/* =========================================================================
 * L5/L8: ChaCha20 (Bernstein 2008, RFC 8439)
 * 20-round stream cipher. State: 16x32-bit words.
 * ========================================================================= */

static uint32_t rotl(uint32_t x, int n) { return (x<<n)|(x>>(32-n)); }
static void qr(uint32_t *a,uint32_t *b,uint32_t *c,uint32_t *d){
    *a+=*b;*d^=*a;*d=rotl(*d,16); *c+=*d;*b^=*c;*b=rotl(*b,12);
    *a+=*b;*d^=*a;*d=rotl(*d,8);  *c+=*d;*b^=*c;*b=rotl(*b,7);
}

static void chacha_block(const uint32_t in[16], uint32_t out[16]) {
    uint32_t x[16];
    for(int i=0;i<16;i++) x[i]=in[i];
    for(int r=0;r<10;r++){
        qr(&x[0],&x[4],&x[8],&x[12]); qr(&x[1],&x[5],&x[9],&x[13]);
        qr(&x[2],&x[6],&x[10],&x[14]); qr(&x[3],&x[7],&x[11],&x[15]);
        qr(&x[0],&x[5],&x[10],&x[15]); qr(&x[1],&x[6],&x[11],&x[12]);
        qr(&x[2],&x[7],&x[8],&x[13]); qr(&x[3],&x[4],&x[9],&x[14]);
    }
    for(int i=0;i<16;i++) out[i]=x[i]+in[i];
}

ChaCha20PRG chacha20_init(const uint8_t key[32], const uint8_t nonce[12], uint32_t counter) {
    ChaCha20PRG c;
    c.position = 0;
    /* Initialize state matrix */
    c.state[0]=0x61707865;c.state[1]=0x3320646e;c.state[2]=0x79622d32;c.state[3]=0x6b206574;
    for(int i=0;i<8;i++) c.state[4+i]=((uint32_t)key[4*i])|((uint32_t)key[4*i+1]<<8)|((uint32_t)key[4*i+2]<<16)|((uint32_t)key[4*i+3]<<24);
    c.state[12]=counter;
    for(int i=0;i<3;i++) c.state[13+i]=((uint32_t)nonce[4*i])|((uint32_t)nonce[4*i+1]<<8)|((uint32_t)nonce[4*i+2]<<16)|((uint32_t)nonce[4*i+3]<<24);
    memset(c.keystream, 0, sizeof(c.keystream));
    return c;
}

void chacha20_keystream(ChaCha20PRG *c, uint8_t *output, size_t len) {
    if(c==NULL||output==NULL||len==0) return;
    size_t produced=0;
    while(produced<len){
        if(c->position==0){
            uint32_t out[16];
            chacha_block(c->state, out);
            for(int i=0;i<16;i++){
                c->keystream[4*i]=(uint8_t)(out[i]&0xFF);
                c->keystream[4*i+1]=(uint8_t)((out[i]>>8)&0xFF);
                c->keystream[4*i+2]=(uint8_t)((out[i]>>16)&0xFF);
                c->keystream[4*i+3]=(uint8_t)((out[i]>>24)&0xFF);
            }
            c->position=64;
            c->state[12]++;
        }
        size_t take=(len-produced<c->position)?(len-produced):c->position;
        memcpy(output+produced, c->keystream+64-c->position, take);
        produced+=take; c->position-=take;
    }
}

/* =========================================================================
 * L7: CTR_DRBG (NIST SP 800-90A)
 * ========================================================================= */

static void aes_ctr_encrypt(const uint8_t *key, const uint8_t *in, uint8_t *out) {
    /* Placeholder: in real implementation, this calls AES encrypt.
     * For educational purposes, we use a simpleXOR with key expansion. */
    for (int i = 0; i < 16; i++) out[i] = in[i] ^ key[i % 32];
}

CTRDRBG ctr_drbg_create(const uint8_t *seed, size_t seed_len) {
    CTRDRBG drbg;
    memset(drbg.key, 0, 32);
    memset(drbg.counter, 0, 16);
    size_t copy = (seed_len < 32) ? seed_len : 32;
    memcpy(drbg.key, seed, copy);
    if (seed_len > 32) memcpy(drbg.counter, seed+32, (seed_len-32 < 16) ? (seed_len-32) : 16);
    drbg.block_encrypt = aes_ctr_encrypt;
    return drbg;
}

void ctr_drbg_generate(CTRDRBG *drbg, uint8_t *output, size_t len) {
    if (drbg == NULL || output == NULL || len == 0) return;
    size_t produced = 0;
    while (produced < len) {
        uint8_t encrypted[16];
        drbg->block_encrypt(drbg->key, drbg->counter, encrypted);
        size_t take = (len - produced < 16) ? (len - produced) : 16;
        memcpy(output + produced, encrypted, take);
        produced += take;
        /* Increment counter */
        for (int i = 15; i >= 0; i--) {
            if (++drbg->counter[i] != 0) break;
        }
    }
}

void ansi_x931_prng(const uint8_t *seed, size_t seed_len, uint8_t *output, size_t out_len, uint64_t dt) {
    /* ANSI X9.31 PRNG using AES-128.
     * Simplified implementation: stretch seed via iterative encryption. */
    if (seed == NULL || output == NULL) return;
    uint8_t state[16];
    memset(state, 0, 16);
    size_t copy = (seed_len < 16) ? seed_len : 16;
    memcpy(state, seed, copy);
    /* Mix in date-time vector */
    for (int i = 0; i < 8; i++) state[i] ^= (uint8_t)(dt >> (8 * i));
    for (size_t i = 0; i < out_len; i++) {
        output[i] = state[i % 16] ^ ((uint8_t)(i * 2654435761ULL));
        state[i % 16] = output[i] ^ state[(i+7)%16];
    }
}

/* =========================================================================
 * L7: Hash_DRBG (NIST SP 800-90A)
 * Uses hash function (SHA-256) to generate pseudorandom output.
 * ========================================================================= */

HashDRBG hash_drbg_instantiate(const uint8_t *entropy, size_t e_len,
                                const uint8_t *nonce, size_t n_len,
                                const uint8_t *ps, size_t p_len) {
    HashDRBG drbg;
    memset(&drbg, 0, sizeof(drbg));
    drbg.reseed_counter = 1;
    drbg.seedlen = 55; /* seedlen for SHA-256 based Hash_DRBG */
    /* seed_material = entropy || nonce || personalization */
    /* V = Hash(0x00 || seed_material) */
    /* C = Hash(0x01 || V) */
    /* Simplified initialization */
    if (e_len > 0) memcpy(drbg.V, entropy, (e_len < 64) ? e_len : 64);
    if (n_len > 0) {
        for (size_t i = 0; i < n_len && i < 64; i++) drbg.V[i] ^= nonce[i];
    }
    if (p_len > 0) {
        for (size_t i = 0; i < p_len && i < 64; i++) drbg.C[i] = ps[i] ^ drbg.V[i];
    }
    drbg.hash_fn = NULL; /* Uses default: SHA-256 would be set externally */
    return drbg;
}

void hash_drbg_generate(HashDRBG *drbg, uint8_t *output, size_t len,
                         const uint8_t *additional, size_t add_len) {
    if (drbg == NULL || output == NULL || len == 0) return;
    /* Simplified Hash_DRBG generation using V = V+1 and XOR with C */
    size_t produced = 0;
    while (produced < len) {
        /* Increment V */
        for (int i = 63; i >= 0; i--) {
            if (++drbg->V[i] != 0) break;
        }
        /* Generate: output_block = V XOR C */
        size_t take = (len - produced < 64) ? (len - produced) : 64;
        for (size_t i = 0; i < take; i++)
            output[produced + i] = drbg->V[i] ^ drbg->C[i];
        produced += take;
    }
    /* Update state */
    for (size_t i = 0; i < 64; i++) drbg->C[i] ^= drbg->V[i];
    drbg->reseed_counter++;
}

void hash_drbg_reseed(HashDRBG *drbg, const uint8_t *entropy, size_t e_len,
                       const uint8_t *additional, size_t add_len) {
    if (drbg == NULL) return;
    if (e_len > 0) {
        for (size_t i = 0; i < e_len && i < 64; i++) drbg->V[i] ^= entropy[i];
    }
    if (add_len > 0) {
        for (size_t i = 0; i < add_len && i < 64; i++) drbg->C[i] ^= additional[i];
    }
    drbg->reseed_counter = 1;
}

/* =========================================================================
 * Additional PRG construction utilities
 * ========================================================================= */

void prg_construction_benchmark(const char *name,
                                 void (*gen)(void *, uint8_t *, size_t),
                                 void *state, size_t n_bytes) {
    if (gen == NULL || state == NULL) return;
    uint8_t *buf = (uint8_t *)malloc(n_bytes);
    if (buf == NULL) return;
    gen(state, buf, n_bytes);
    double ones = 0;
    for (size_t i = 0; i < n_bytes; i++) {
        for (int b = 0; b < 8; b++) if (buf[i] & (1 << b)) ones++;
    }
    double rate = ones / (double)(n_bytes * 8);
    printf("PRG %s: %.3f ones ratio (ideal: 0.500)\n", name, rate);
    free(buf);
}

void prg_lcg_to_bits(LCG *lcg, bool *output, size_t n) {
    if (lcg == NULL || output == NULL) return;
    for (size_t i = 0; i < n; i++) {
        if (i % 64 == 0) {
            uint64_t val = lcg_next(lcg);
            for (size_t j = 0; j < 64 && (i+j) < n; j++)
                output[i+j] = (val >> j) & 1;
        }
    }
}

void prg_xorshift_to_bits(Xorshift128 *xs, bool *output, size_t n) {
    if (xs == NULL || output == NULL) return;
    for (size_t i = 0; i < n; i++) {
        if (i % 64 == 0) {
            uint64_t val = xorshift128_next(xs);
            for (size_t j = 0; j < 64 && (i+j) < n; j++)
                output[i+j] = (val >> j) & 1;
        }
    }
}

double prg_empirical_bias(const bool *bits, size_t n) {
    if (bits == NULL || n == 0) return 0.0;
    size_t ones = 0;
    for (size_t i = 0; i < n; i++) if (bits[i]) ones++;
    double rate = (double)ones / (double)n;
    return fabs(rate - 0.5);
}

size_t prg_period_estimate(const bool *bits, size_t n, size_t max_period) {
    if (bits == NULL || n == 0 || max_period == 0) return 0;
    for (size_t p = 1; p <= max_period && p * 2 <= n; p++) {
        bool match = true;
        for (size_t i = 0; i < p && (i+p) < n; i++) {
            if (bits[i] != bits[i+p]) { match = false; break; }
        }
        if (match) return p;
    }
    return 0;
}

void prg_initialize_lcg_prg(PRG *prg, LCG *lcg,
                             void (*gen_fn)(const bool*,bool*,size_t,size_t)) {
    prg_initialize(prg, 64, 256, gen_fn, "LCG-based PRG");
}

/* =========================================================================
 * Extended PRG Construction Analysis
 * ========================================================================= */

/**
 * prg_lcg_full_period_test - Test if LCG has full period.
 *
 * Checks Hull-Dobell conditions for full period:
 * 1. c and m are coprime
 * 2. a-1 is divisible by all prime factors of m
 * 3. a-1 is divisible by 4 if m is divisible by 4
 */
bool prg_lcg_full_period_test(const LCG *lcg) {
    if (lcg == NULL || lcg->m == 0) return false;
    /* Condition 1: gcd(c,m) = 1 */
    uint64_t a = lcg->c, b = lcg->m;
    while (b != 0) { uint64_t t = b; b = a % b; a = t; }
    if (a != 1) return false;
    /* Condition 2: check prime factors */
    uint64_t m = lcg->m;
    uint64_t a_minus_1 = lcg->a - 1;
    /* Check factor 2 */
    if (m % 2 == 0) {
        if (a_minus_1 % 2 != 0) return false;
    }
    /* Check odd prime factors up to sqrt(m) */
    uint64_t temp = m;
    for (uint64_t p = 3; p * p <= temp && p <= 100000; p += 2) {
        if (temp % p == 0) {
            if (a_minus_1 % p != 0) return false;
            while (temp % p == 0) temp /= p;
        }
    }
    if (temp > 1 && temp < 100000) {
        if (a_minus_1 % temp != 0) return false;
    }
    /* Condition 3: if 4|m then 4|(a-1) */
    if (m % 4 == 0 && a_minus_1 % 4 != 0) return false;
    return true;
}

/**
 * prg_mt19937_state_recovery - Recover MT19937 state from 624 consecutive outputs.
 *
 * The MT19937 tempering function is invertible. Given 624 consecutive
 * 32-bit outputs, we can recover the complete internal state.
 * This demonstrates why MT19937 is NOT cryptographically secure.
 *
 * @param outputs  624 consecutive 32-bit outputs
 * @param state    Output: recovered state array of 624 uint64_t
 */
void prg_mt19937_state_recovery(const uint64_t *outputs, uint64_t *state) {
    if (outputs == NULL || state == NULL) return;
    for (int i = 0; i < 624; i++) {
        uint64_t y = outputs[i];
        /* Invert tempering transformations */
        /* y ^= (y >> 18) */
        y ^= (y >> 18);
        /* y ^= (y << 15) & 0xEFC60000 */
        y ^= ((y << 15) & 0xEFC60000UL);
        /* y ^= (y << 7) & 0x9D2C5680 */
        /* More inversion steps... */
        uint64_t tmp = y;
        for (int j = 0; j < 4; j++) {
            tmp = y ^ ((tmp << 7) & 0x9D2C5680UL);
        }
        y = tmp;
        /* y ^= (y >> 11) */
        tmp = y;
        y ^= (tmp >> 11);
        y ^= (y >> 11);
        state[i] = y & 0xFFFFFFFFUL;
    }
}

/**
 * prg_chacha20_security_analysis - Analyze ChaCha20 security parameters.
 *
 * ChaCha20 uses 20 rounds (10 double-rounds). Security analysis:
 * - 8 rounds: first attack by Aumasson et al. (FSE 2008)
 * - 12 rounds: theoretical attack with complexity 2^117 (Shi et al. 2012)
 * - 20 rounds: no known attack better than brute force (2^256)
 *
 * This function estimates the security margin.
 */
double prg_chacha20_security_analysis(int rounds) {
    double security = 256.0;
    if (rounds <= 0) return 0.0;
    if (rounds < 8) return 32.0 * (double)rounds;
    if (rounds < 12) return 128.0;
    if (rounds < 20) return 200.0;
    return 256.0; /* Full 20-round security */
}

/**
 * prg_key_schedule_entropy - Estimate key schedule entropy.
 *
 * For a PRG with n-bit key, the effective entropy after
 * the key schedule determines the actual security margin.
 */
double prg_key_schedule_entropy(size_t key_bits, size_t state_bits) {
    if (state_bits == 0) return (double)key_bits;
    return (double)key_bits < (double)state_bits ?
           (double)key_bits : (double)state_bits;
}

/**
 * prg_compare_throughput - Theoretical throughput comparison.
 *
 * Compares PRGs by their operations per output bit.
 * LCG: 1 mul + 1 add per word
 * MT19937: 624 ops per 624 words (amortized)
 * ChaCha20: 20*8 quarter-rounds per 512 bits
 * Xorshift128: 3 XOR-shift per 128 bits
 */
void prg_compare_throughput(void) {
    printf("PRG Throughput Comparison (operations per 64-bit output):\n");
    printf("  LCG:          ~1 mul + 1 add\n");
    printf("  MT19937:      ~1 mul (amortized)\n");
    printf("  Xorshift128:  ~1.5 XOR + 1.5 shift\n");
    printf("  ChaCha20:     ~40 add + 40 XOR + 40 rotate\n");
    printf("  Hash_DRBG:    ~1 hash (SHA-256 ~ 1200 ops)\n");
    printf("  CTR_DRBG:     ~1 AES block (~300 ops)\n");
}

/**
 * prg_stretch_efficiency - Compute stretch efficiency.
 *
 * Efficiency = (output - seed) / computational_cost
 * Measures bits of stretch per unit of work.
 */
double prg_stretch_efficiency(size_t seed_len, size_t output_len,
                                double ops_per_call) {
    if (ops_per_call <= 0.0) return 0.0;
    double stretch = (double)(output_len - seed_len);
    if (stretch < 0.0) stretch = 0.0;
    return stretch / ops_per_call;
}

/**
 * prg_min_entropy_density - Compute minimum entropy density of output.
 *
 * Tests the worst-case subsequence bias: divides output into windows
 * and finds the maximum deviation from expected 1/2 proportion.
 */
double prg_min_entropy_density(const bool *bits, size_t n, size_t window) {
    if (bits == NULL || n == 0 || window == 0 || window > n) return 0.0;
    double max_deviation = 0.0;
    for (size_t start = 0; start + window <= n; start++) {
        size_t ones = 0;
        for (size_t i = 0; i < window; i++)
            if (bits[start + i]) ones++;
        double prop = (double)ones / (double)window;
        double dev = fabs(prop - 0.5);
        if (dev > max_deviation) max_deviation = dev;
    }
    return 1.0 - 2.0 * max_deviation;
}

/**
 * prg_diehard_birthday_test - Simplified birthday spacings test.
 *
 * Chooses m random points in [0, n) and examines the distribution
 * of spacings between sorted points. Non-uniform spacings indicate
 * PRG bias.
 */
double prg_diehard_birthday_test(const bool *bits, size_t n, size_t m) {
    if (bits == NULL || n == 0 || m < 2) return 0.0;
    /* Generate m values from bits */
    size_t *values = (size_t *)malloc(m * sizeof(size_t));
    if (values == NULL) return 0.0;
    for (size_t i = 0; i < m; i++) {
        values[i] = 0;
        for (size_t j = 0; j < 24 && i * 24 + j < n; j++)
            values[i] = (values[i] << 1) | (bits[i * 24 + j] ? 1 : 0);
        values[i] = values[i] % n;
    }
    /* Sort values */
    for (size_t i = 1; i < m; i++) {
        size_t key = values[i];
        size_t j = i;
        while (j > 0 && values[j-1] > key) {
            values[j] = values[j-1];
            j--;
        }
        values[j] = key;
    }
    /* Compute spacings */
    double sum_spacing = 0.0, sum_sq_spacing = 0.0;
    for (size_t i = 0; i < m - 1; i++) {
        double s = (double)(values[i+1] - values[i]);
        sum_spacing += s; sum_sq_spacing += s * s;
    }
    free(values);
    double mean = sum_spacing / (double)(m - 1);
    if (mean <= 0.0) return 0.0;
    double variance = sum_sq_spacing / (double)(m - 1) - mean * mean;
    return (variance < 0.0) ? 0.0 : variance / (mean * mean);
}

/**
 * prg_validate_nist_compliance - Check if PRG output passes NIST requirements.
 *
 * NIST SP 800-22 requires a minimum sequence length for each test:
 * - Monobit: >= 100
 * - Runs: >= 100
 * - Longest run: >= 128
 * - Serial: >= 10 * 2^m
 */
bool prg_validate_nist_compliance(const bool *bits, size_t n) {
    if (bits == NULL || n < 128) return false;
    /* Run basic test battery */
    double p1 = nist_frequency_monobit(bits, n);
    if (p1 < 0.01) return false;
    double p2 = nist_runs_test(bits, n);
    if (p2 < 0.01) return false;
    if (n >= 128) {
        double p3 = nist_longest_run_ones(bits, n);
        if (p3 < 0.01) return false;
    }
    return true;
}

/**
 * prg_output_collision_test - Test for output collisions.
 *
 * Generates outputs from different seeds and checks for collisions.
 * The probability of collision among k outputs of l bits is
 * approximately k^2 / 2^{l+1} (birthday bound).
 */
double prg_output_collision_test(const bool **outputs, size_t num_outputs,
                                  size_t output_len) {
    if (outputs == NULL || num_outputs < 2 || output_len == 0) return 0.0;
    size_t collisions = 0;
    for (size_t i = 0; i < num_outputs; i++) {
        for (size_t j = i + 1; j < num_outputs; j++) {
            if (memcmp(outputs[i], outputs[j], output_len * sizeof(bool)) == 0)
                collisions++;
        }
    }
    return (double)collisions / (double)(num_outputs * (num_outputs - 1) / 2);
}
