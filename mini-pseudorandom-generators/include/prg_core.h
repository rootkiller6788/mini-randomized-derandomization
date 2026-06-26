#ifndef PRG_CORE_H
#define PRG_CORE_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>

typedef struct { size_t seed_len, output_len, stretch; void (*gen)(const bool*,bool*,size_t,size_t); double adv; const char *name; } PRG;
typedef struct { bool (*pred)(const bool*,size_t); size_t queries; double adv; } Distinguisher;
typedef struct { bool (*pred_next)(const bool*,size_t); double prob; size_t pos; } NextBitPredictor;

/* L2: Statistical distance */
double prg_statistical_distance(const double *p, const double *q, size_t n);
double prg_distinguish_advantage(const bool *g_out, const bool *u_out, size_t n, size_t trials);
bool prg_yao_theorem_verify(const bool *prefix, size_t prefix_len, bool predicted_bit, bool actual_bit, size_t trials);
size_t prg_bpp_trials_for_error(double delta);
bool prg_has_stretch(const PRG *prg);
double prg_stretch_ratio(const PRG *prg);
size_t prg_cycle_detect(const PRG *prg, const bool *seed, size_t *period);
double prg_hybrid_argument_bound(size_t m, double nbp_adv);
void prg_stream_encrypt(const bool *key, size_t kl, const bool *pt, size_t ml, bool *ct, const PRG *prg);
bool prg_derandomize_bpp(bool (*dec)(const char*,size_t,const bool*,size_t), const char *inp, size_t ilen, const PRG *prg);

/* Extended L2-L7 functions */
double prg_computational_distance(const bool *seq_a, const bool *seq_b, size_t n);
void prg_seed_expand(const PRG *prg, const bool *seed, bool *output);
void prg_compose_stretch(const PRG *prg, const bool *seed, size_t k, bool *output);
double prg_iw97_bound(size_t m, double epsilon, size_t input_len);
double prg_security_level(size_t seed_len, double advantage);
size_t prg_min_seed_length(size_t target_security, double advantage);
bool prg_forward_secrecy_check(const bool *prev_output, const bool *curr_state, size_t seed_len);
bool prg_is_negligible(double advantage, size_t n);
bool prg_compose(PRG *prgout, const PRG *g1, const PRG *g2);
size_t prg_bit_generation_rate(const PRG *prg);
size_t prg_entropy_loss(const PRG *prg);
void prg_initialize(PRG *prg, size_t seed_len, size_t output_len, void (*gen)(const bool*,bool*,size_t,size_t), const char *name);
void prg_create_distinguisher(Distinguisher *d, bool (*pred)(const bool*,size_t), double adv);
void prg_create_predictor(NextBitPredictor *nbp, bool (*pred)(const bool*,size_t), double prob, size_t pos);
void prg_print_info(const PRG *prg);

/* Nisan-Wigderson framework */
size_t prg_nw_design_size(size_t m, size_t k, size_t r);
size_t prg_nw_output_length(size_t d, size_t r);
size_t prg_nw_hardness_to_stretch(size_t hardness, size_t seed_length);
double prg_security_parameter(size_t seed_len, size_t time_bound);
double prg_advantage_decay(double initial_adv, size_t trials);
double prg_unpredictability_to_indistinguishability(double nbp_advantage, size_t num_bits);
double prg_indistinguishability_to_unpredictability(double distinguish_adv, size_t num_bits);
void prg_quantitative_yao(double nbp_adv, double dist_adv, size_t num_bits, double *out_nbp_bound, double *out_dist_bound);

/* Extended PRG analysis */
double prg_entropy_estimate(const bool *bits, size_t n, size_t block_size);
double prg_compression_ratio(const bool *bits, size_t n);
double prg_correlation_matrix(const bool *bits, size_t n, size_t num_positions);
double prg_mutual_information(const bool *bits, size_t n);
double prg_chi_squared_uniformity(const bool *bits, size_t n, size_t buckets);
double prg_poker_test(const bool *bits, size_t n, size_t m);

#endif /* PRG_CORE_H */
