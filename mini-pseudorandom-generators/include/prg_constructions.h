#ifndef PRG_CONSTRUCTIONS_H
#define PRG_CONSTRUCTIONS_H
#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include "prg_core.h"

typedef struct { uint64_t a, c, m, state; } LCG;
LCG lcg_create(uint64_t seed, uint64_t a, uint64_t c, uint64_t m);
uint64_t lcg_next(LCG *lcg);
double lcg_quality_test(const LCG *lcg, size_t n);

typedef struct {
    uint64_t state[624];
    int index;
    uint64_t seed;
} MT19937;
MT19937 mt19937_create(uint64_t seed);
uint64_t mt19937_extract(MT19937 *mt);
void mt19937_generate_bits(MT19937 *mt, bool *output, size_t n);

typedef struct { uint64_t s[2]; } Xorshift128;
Xorshift128 xorshift128_create(uint64_t seed1, uint64_t seed2);
uint64_t xorshift128_next(Xorshift128 *xs);

typedef struct {
    uint32_t state[16];
    uint32_t keystream[16];
    size_t position;
} ChaCha20PRG;
ChaCha20PRG chacha20_init(const uint8_t key[32], const uint8_t nonce[12], uint32_t counter);
void chacha20_keystream(ChaCha20PRG *c, uint8_t *output, size_t len);

typedef struct {
    uint8_t key[32];
    uint8_t counter[16];
    void (*block_encrypt)(const uint8_t *key, const uint8_t *in, uint8_t *out);
} CTRDRBG;
CTRDRBG ctr_drbg_create(const uint8_t *seed, size_t seed_len);
void ctr_drbg_generate(CTRDRBG *drbg, uint8_t *output, size_t len);
void ansi_x931_prng(const uint8_t *seed, size_t seed_len, uint8_t *output, size_t out_len, uint64_t dt);

typedef struct {
    uint8_t V[64], C[64];
    uint64_t reseed_counter;
    void (*hash_fn)(const uint8_t *in, size_t len, uint8_t *out);
    size_t seedlen;
} HashDRBG;
HashDRBG hash_drbg_instantiate(const uint8_t *entropy, size_t e_len, const uint8_t *nonce, size_t n_len, const uint8_t *personalization, size_t p_len);
void hash_drbg_generate(HashDRBG *drbg, uint8_t *output, size_t len, const uint8_t *additional, size_t add_len);
void hash_drbg_reseed(HashDRBG *drbg, const uint8_t *entropy, size_t e_len, const uint8_t *additional, size_t add_len);

/* Extended PRG construction analysis */
bool prg_lcg_full_period_test(const LCG *lcg);
void prg_mt19937_state_recovery(const uint64_t *outputs, uint64_t *state);
double prg_chacha20_security_analysis(int rounds);
double prg_key_schedule_entropy(size_t key_bits, size_t state_bits);
void prg_compare_throughput(void);
double prg_stretch_efficiency(size_t seed_len, size_t output_len, double ops_per_call);
double prg_min_entropy_density(const bool *bits, size_t n, size_t window);
double prg_diehard_birthday_test(const bool *bits, size_t n, size_t m);
bool prg_validate_nist_compliance(const bool *bits, size_t n);
double prg_output_collision_test(const bool **outputs, size_t num_outputs, size_t output_len);

/* Utility functions */
void prg_construction_benchmark(const char *name, void (*gen)(void*,uint8_t*,size_t), void *state, size_t n_bytes);
void prg_lcg_to_bits(LCG *lcg, bool *output, size_t n);
void prg_xorshift_to_bits(Xorshift128 *xs, bool *output, size_t n);
double prg_empirical_bias(const bool *bits, size_t n);
size_t prg_period_estimate(const bool *bits, size_t n, size_t max_period);
void prg_initialize_lcg_prg(PRG *prg, LCG *lcg, void (*gen_fn)(const bool*,bool*,size_t,size_t));
#endif
