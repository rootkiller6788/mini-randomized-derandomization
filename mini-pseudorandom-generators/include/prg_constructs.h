#ifndef PRG_CONSTRUCTS_H
#define PRG_CONSTRUCTS_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>

typedef struct { uint64_t a, c, m, state; } LCGState;
typedef struct { uint64_t s[624]; int idx; } MTState;
typedef struct { uint64_t s[2]; } XS128State;
typedef struct { uint32_t s[16]; uint32_t ks[16]; size_t pos; } ChaChaState;
typedef struct { uint8_t V[64], C[64]; uint64_t ctr; } HashDRBGState;

LCGState lcg_init(uint64_t seed, uint64_t a, uint64_t c, uint64_t m);
uint64_t lcg_next(LCGState *s);
double lcg_quality(const LCGState *s, size_t n);
MTState mt19937_init(uint64_t seed);
uint64_t mt19937_next(MTState *s);
void mt19937_fill(MTState *s, bool *out, size_t n);
XS128State xs128_init(uint64_t s0, uint64_t s1);
uint64_t xs128_next(XS128State *s);
ChaChaState chacha20_init(const uint8_t key[32], const uint8_t nonce[12], uint32_t ctr);
void chacha20_keystream_bytes(ChaChaState *s, uint8_t *out, size_t n);
HashDRBGState hash_drbg_init(const uint8_t *seed, size_t slen);
void hash_drbg_gen(HashDRBGState *s, uint8_t *out, size_t n);
void hash_drbg_reseed(HashDRBGState *s, const uint8_t *seed, size_t slen);
void ansi_x931_prng(const uint8_t *seed, size_t slen, uint8_t *out, size_t olen, uint64_t dt);
void ctr_drbg_gen(const uint8_t *key, const uint8_t *ctr, uint8_t *out, size_t n);

#endif /* PRG_CONSTRUCTS_H */
