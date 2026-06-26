#ifndef HARDCORE_PREDICATES_H
#define HARDCORE_PREDICATES_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>

typedef struct { bool (*pred)(const bool*,size_t); bool (*owf)(const bool*,size_t,bool*,size_t); double hardness; } HardcorePred;
typedef struct { bool (*eval)(const bool*,size_t,bool*,size_t); size_t ilen, olen; double inv_prob; } OneWayFn;
bool hc_goldreich_levin(const bool *x, const bool *r, size_t n, bool *result);
bool hc_blum_micali_msb(uint64_t x, uint64_t g, uint64_t p, bool *bit);
bool hc_blum_blum_shub(uint64_t seed, uint64_t N, size_t count, bool *bits);
double hc_advantage_bound(size_t n, size_t queries, double epsilon);
bool hc_verify(const HardcorePred *hcp, size_t trials, double *empirical_adv);
HardcorePred hc_construct_gl(const OneWayFn *owf);
bool owf_test_inversion(const OneWayFn *owf, size_t trials, double *rate);
OneWayFn owf_rsa_create(uint64_t n, uint64_t e);
OneWayFn owf_dlog_create(uint64_t p, uint64_t g);
uint64_t modpow_u64(uint64_t base, uint64_t exp, uint64_t mod);
bool is_prime_u64(uint64_t n);

/* Extended hardcore predicate analysis */
bool hc_combine(const bool *bits, size_t k, bool *result);
bool hc_universal_one_way(const bool *x, size_t n, const bool *r, size_t k, bool *hash);
double hc_statistical_test(const bool *x, size_t n, size_t samples);
double hc_amplify(double epsilon, size_t k);
size_t hc_security_to_keylen(size_t lambda, const char *scheme);
size_t hc_list_decode_gl(size_t n, bool (*oracle)(const bool*,size_t), bool **list, size_t max_list);
double hc_yao_xor_lemma(double epsilon, size_t k);
size_t hc_levin_one_way(const char *description, const char *input, size_t time_bound);
double hc_impagliazzo_hardcore(size_t n, double delta, size_t S);
void hc_generic_construction(const bool *x, size_t n, bool *output, void (*f)(const bool*,size_t,bool*,size_t), bool (*b)(const bool*,size_t));
double hc_quantitative_security(size_t n, double owf_advantage);
void hc_hybrid_prg_construction(const bool *x, size_t n, bool *output, size_t k, bool (**preds)(const bool*,size_t));
double hc_small_bias_test(bool (*pred)(const bool*,size_t), const bool **inputs, size_t num_inputs, size_t n);
bool hc_unpredictable_sampling(bool (*adversary)(void), size_t trials, double threshold);

#endif /* HARDCORE_PREDICATES_H */
