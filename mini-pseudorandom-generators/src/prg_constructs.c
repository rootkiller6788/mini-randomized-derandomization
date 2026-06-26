#include "prg_constructs.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* =========================================================================
 * L3: Lightweight PRG State Implementations
 *
 * This file provides simpler, education-oriented implementations
 * of PRG state machines. These complement prg_constructions.c which
 * provides full industrial-grade implementations.
 *
 * Focus: Mathematical structure, state transitions, and analysis.
 * ========================================================================= */

/* LCGState - Minimal Linear Congruential Generator */

LCGState lcg_init(uint64_t seed, uint64_t a, uint64_t c, uint64_t m) {
    LCGState s;
    s.a = a; s.c = c; s.m = m;
    s.state = (m > 0) ? seed % m : seed;
    if (s.state == 0 && c == 0) s.state = 1;
    return s;
}

uint64_t lcg_next(LCGState *s) {
    if (s == NULL) return 0;
    if (s->m > 0 && s->m <= 0xFFFFFFFFULL) {
        s->state = ((s->a % s->m) * (s->state % s->m) + s->c) % s->m;
    } else {
        s->state = s->a * s->state + s->c;
        if (s->m > 0) s->state %= s->m;
    }
    return s->state;
}

double lcg_quality(const LCGState *s, size_t n) {
    if (s == NULL || n == 0) return 0.0;
    LCGState copy = *s;
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) sum += (double)lcg_next(&copy);
    double mean = sum / (double)n;
    double expected = (s->m == 0) ? 9.22e18 : (double)s->m / 2.0;
    double dev = expected > 0 ? fabs(mean - expected) / expected : 1.0;
    return (1.0 - dev) < 0.0 ? 0.0 : (1.0 - dev);
}

/* MTState - Lightweight Mersenne Twister wrapper */

#define MT_N 624
#define MT_M 397

MTState mt19937_init(uint64_t seed) {
    MTState mt; mt.idx = MT_N;
    mt.s[0] = seed & 0xFFFFFFFFUL;
    for (int i = 1; i < MT_N; i++) {
        mt.s[i] = (1812433253UL * (mt.s[i-1] ^ (mt.s[i-1] >> 30)) + i);
        mt.s[i] &= 0xFFFFFFFFUL;
    }
    return mt;
}

static void mt_twist(MTState *mt) {
    static const uint64_t mag[2] = {0UL, 0x9908B0DFUL};
    for (int i = 0; i < MT_N; i++) {
        uint64_t x = (mt->s[i] & 0x80000000UL) | (mt->s[(i+1)%MT_N] & 0x7FFFFFFFUL);
        mt->s[i] = mt->s[(i+MT_M)%MT_N] ^ (x>>1) ^ mag[x&1UL];
    }
    mt->idx = 0;
}

uint64_t mt19937_next(MTState *s) {
    if (s == NULL) return 0;
    if (s->idx >= MT_N) mt_twist(s);
    uint64_t y = s->s[s->idx++];
    y ^= (y>>11); y ^= (y<<7) & 0x9D2C5680UL;
    y ^= (y<<15) & 0xEFC60000UL; y ^= (y>>18);
    return y & 0xFFFFFFFFUL;
}

void mt19937_fill(MTState *s, bool *out, size_t n) {
    if (s == NULL || out == NULL) return;
    for (size_t i = 0; i < n; i++) {
        if (i % 32 == 0) {
            uint64_t v = mt19937_next(s);
            for (size_t j = 0; j < 32 && (i+j) < n; j++)
                out[i+j] = (v >> j) & 1;
        }
    }
}

/* XS128State - Xorshift128 lightweight implementation */

XS128State xs128_init(uint64_t s0, uint64_t s1) {
    XS128State xs;
    xs.s[0] = (s0 == 0) ? 123456789ULL : s0;
    xs.s[1] = (s1 == 0) ? 987654321ULL : s1;
    return xs;
}

uint64_t xs128_next(XS128State *s) {
    if (s == NULL) return 0;
    uint64_t x = s->s[0], y = s->s[1];
    uint64_t result = x + y;
    s->s[0] = y;
    x ^= x << 23; x ^= x >> 17;
    x ^= y; x ^= y >> 26;
    s->s[1] = x;
    return result;
}

/* ChaChaState - Minimal ChaCha20 state machine */

static uint32_t rotl32(uint32_t v, int n) { return (v<<n)|(v>>(32-n)); }

static void chacha_qr(uint32_t *a, uint32_t *b, uint32_t *c, uint32_t *d) {
    *a+=*b;*d^=*a;*d=rotl32(*d,16);*c+=*d;*b^=*c;*b=rotl32(*b,12);
    *a+=*b;*d^=*a;*d=rotl32(*d,8);*c+=*d;*b^=*c;*b=rotl32(*b,7);
}

ChaChaState chacha20_init(const uint8_t key[32], const uint8_t nonce[12], uint32_t ctr) {
    ChaChaState cs;
    cs.pos = 0;
    cs.s[0]=0x61707865;cs.s[1]=0x3320646e;cs.s[2]=0x79622d32;cs.s[3]=0x6b206574;
    for(int i=0;i<8;i++) cs.s[4+i]=((uint32_t)key[4*i])|((uint32_t)key[4*i+1]<<8)|((uint32_t)key[4*i+2]<<16)|((uint32_t)key[4*i+3]<<24);
    cs.s[12]=ctr;
    for(int i=0;i<3;i++) cs.s[13+i]=((uint32_t)nonce[4*i])|((uint32_t)nonce[4*i+1]<<8)|((uint32_t)nonce[4*i+2]<<16)|((uint32_t)nonce[4*i+3]<<24);
    memset(cs.ks,0,sizeof(cs.ks));
    return cs;
}

static void chacha20_block_core(const uint32_t in[16], uint32_t out[16]) {
    uint32_t x[16];
    for(int i=0;i<16;i++)x[i]=in[i];
    for(int r=0;r<10;r++){
        chacha_qr(&x[0],&x[4],&x[8],&x[12]);chacha_qr(&x[1],&x[5],&x[9],&x[13]);
        chacha_qr(&x[2],&x[6],&x[10],&x[14]);chacha_qr(&x[3],&x[7],&x[11],&x[15]);
        chacha_qr(&x[0],&x[5],&x[10],&x[15]);chacha_qr(&x[1],&x[6],&x[11],&x[12]);
        chacha_qr(&x[2],&x[7],&x[8],&x[13]);chacha_qr(&x[3],&x[4],&x[9],&x[14]);
    }
    for(int i=0;i<16;i++)out[i]=x[i]+in[i];
}

void chacha20_keystream_bytes(ChaChaState *s, uint8_t *out, size_t n) {
    if(s==NULL||out==NULL||n==0)return;
    size_t produced=0;
    while(produced<n){
        if(s->pos==0){
            uint32_t block[16];
            chacha20_block_core(s->s,block);
            for(int i=0;i<16;i++){
                s->ks[4*i]=(uint8_t)block[i];
                s->ks[4*i+1]=(uint8_t)(block[i]>>8);
                s->ks[4*i+2]=(uint8_t)(block[i]>>16);
                s->ks[4*i+3]=(uint8_t)(block[i]>>24);
            }
            s->pos=64;s->s[12]++;
        }
        size_t take=(n-produced<(64-s->pos))?(n-produced):(64-s->pos);
        memcpy(out+produced,s->ks+64-s->pos,take);
        produced+=take;s->pos-=take;
    }
}

/* HashDRBGState - Simple Hash-based DRBG state machine */

HashDRBGState hash_drbg_init(const uint8_t *seed, size_t slen) {
    HashDRBGState drbg;
    memset(&drbg,0,sizeof(drbg));
    drbg.ctr=1;
    if(slen>0)memcpy(drbg.V,seed,(slen<64)?slen:64);
    for(size_t i=0;i<64;i++)drbg.C[i]=drbg.V[(i+13)%64]^0xA5;
    return drbg;
}

void hash_drbg_gen(HashDRBGState *s, uint8_t *out, size_t n) {
    if(s==NULL||out==NULL||n==0)return;
    size_t produced=0;
    while(produced<n){
        for(int i=63;i>=0;i--){if(++s->V[i]!=0)break;}
        size_t take=(n-produced<64)?(n-produced):64;
        for(size_t i=0;i<take;i++)out[produced+i]=s->V[i]^s->C[i];
        produced+=take;
    }
    for(size_t i=0;i<64;i++)s->C[i]^=s->V[i];
    s->ctr++;
}

void hash_drbg_reseed(HashDRBGState *s, const uint8_t *seed, size_t slen) {
    if(s==NULL)return;
    if(slen>0){for(size_t i=0;i<slen&&i<64;i++)s->V[i]^=seed[i];}
    s->ctr=1;
}

/* ANSI X9.31 PRNG */
void ansi_x931_prng(const uint8_t *seed, size_t slen, uint8_t *out, size_t olen, uint64_t dt) {
    if(seed==NULL||out==NULL)return;
    uint8_t state[16];
    memset(state,0,16);
    size_t cp=(slen<16)?slen:16;
    memcpy(state,seed,cp);
    for(int i=0;i<8;i++)state[i]^=(uint8_t)(dt>>(8*i));
    for(size_t i=0;i<olen;i++){
        out[i]=state[i%16]^((uint8_t)(i*2654435761ULL));
        state[i%16]=out[i]^state[(i+7)%16];
    }
}

/* CTR_DRBG simple implementation */
void ctr_drbg_gen(const uint8_t *key, const uint8_t *ctr, uint8_t *out, size_t n) {
    if(key==NULL||ctr==NULL||out==NULL||n==0)return;
    uint8_t counter[16];
    memcpy(counter,ctr,16);
    for(size_t i=0;i<n;i++){
        out[i]=key[i%32]^counter[i%16]^((uint8_t)i);
        for(int j=15;j>=0;j--){if(++counter[j]!=0)break;}
    }
}

/* =========================================================================
 * Additional analysis functions for PRG state machines
 * ========================================================================= */

/**
 * lcg_period_estimate - Estimate LCG period via cycle detection.
 */
size_t lcg_period_estimate(LCGState *s, size_t max_steps) {
    if (s == NULL || max_steps == 0) return 0;
    LCGState tortoise = *s;
    LCGState hare = *s;
    for (size_t step = 0; step < max_steps; step++) {
        lcg_next(&tortoise);
        lcg_next(&hare); lcg_next(&hare);
        if (tortoise.state == hare.state) {
            /* Find period */
            size_t period = 1;
            LCGState runner = tortoise;
            lcg_next(&runner);
            while (runner.state != tortoise.state && period < max_steps) {
                lcg_next(&runner);
                period++;
            }
            return period;
        }
    }
    return 0;
}

/**
 * xs128_period_confirm - Confirm xorshift128 period (should be 2^128-1).
 */
bool xs128_period_confirm(XS128State *s, size_t max_steps) {
    if (s == NULL || max_steps == 0) return false;
    XS128State initial = *s;
    for (size_t step = 1; step < max_steps; step++) {
        xs128_next(s);
        if (s->s[0] == initial.s[0] && s->s[1] == initial.s[1])
            return true;
    }
    return false;
}

/**
 * prg_state_entropy - Estimate entropy of PRG state in bits.
 */
double prg_state_entropy(const uint64_t *state, size_t words) {
    if (state == NULL || words == 0) return 0.0;
    /* Count unique bit patterns approximately via popcount */
    size_t total_bits = 0;
    for (size_t i = 0; i < words; i++) {
        uint64_t v = state[i];
        while (v) { total_bits += (v & 1); v >>= 1; }
    }
    double p = (double)total_bits / (double)(words * 64);
    if (p <= 0.0 || p >= 1.0) return 0.0;
    return -p * log2(p) - (1.0-p) * log2(1.0-p);
}

/**
 * prg_state_distance - Hamming distance between two PRG states.
 */
size_t prg_state_distance(const uint64_t *a, const uint64_t *b, size_t words) {
    if (a == NULL || b == NULL || words == 0) return 0;
    size_t dist = 0;
    for (size_t i = 0; i < words; i++) {
        uint64_t diff = a[i] ^ b[i];
        while (diff) { dist += (diff & 1); diff >>= 1; }
    }
    return dist;
}

/**
 * prg_validate_seed - Check if a seed value is appropriate for PRG use.
 *
 * Returns true if seed has at least some entropy (not all zeros/ones).
 */
bool prg_validate_seed(const uint8_t *seed, size_t len) {
    if (seed == NULL || len == 0) return false;
    bool all_same = true;
    uint8_t first_byte = seed[0];
    for (size_t i = 1; i < len; i++) {
        if (seed[i] != first_byte) { all_same = false; break; }
    }
    /* Reject all-zeros or all-ones seeds */
    if (all_same && (first_byte == 0 || first_byte == 0xFF)) return false;
    /* Reject zero-entropy constant patterns */
    size_t transitions = 0;
    for (size_t i = 1; i < len; i++) {
        if (seed[i] != seed[i-1]) transitions++;
    }
    /* At least some bit changes */
    if (transitions < len / 16 && len > 16) return false;
    return true;
}

/**
 * prg_mix_seed - Mix entropy sources into a seed.
 *
 * Combines multiple entropy sources using XOR-then-shift.
 */
void prg_mix_seed(uint8_t *seed, size_t len,
                   const uint8_t *source1, size_t len1,
                   const uint8_t *source2, size_t len2) {
    if (seed == NULL || len == 0) return;
    memset(seed, 0, len);
    if (source1 != NULL) {
        for (size_t i = 0; i < len1 && i < len; i++)
            seed[i % len] ^= source1[i];
    }
    if (source2 != NULL) {
        for (size_t i = 0; i < len2 && i < len; i++)
            seed[(i + 7) % len] ^= source2[i];
    }
    /* Diffusion: mix bytes */
    for (size_t i = 1; i < len; i++)
        seed[i] ^= seed[i-1] * 181 + 37;
}

/**
 * prg_generate_from_state - Generic PRG bit generation from state.
 */
void prg_generate_from_state(const uint8_t *state, size_t state_len,
                              bool *output, size_t out_len) {
    if (state == NULL || output == NULL || state_len == 0 || out_len == 0) return;
    for (size_t i = 0; i < out_len; i++) {
        output[i] = (state[i % state_len] >> (i % 8)) & 1;
    }
}
