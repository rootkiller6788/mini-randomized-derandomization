/**
 * space_derand.h - Core Space Derandomization: Definitions and Structures
 */

#ifndef SPACE_DERAND_H
#define SPACE_DERAND_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>

typedef struct {
    size_t space;
    size_t time;
    size_t num_states;
    size_t alphabet_size;
    int    initial_state;
    size_t num_accept;
    int   *accept_states;
    void  *transition_data;
    size_t trans_size;
} SpaceTM;

typedef struct {
    int    state;
    size_t input_head;
    size_t work_head;
    char  *work_tape;
    size_t space_used;
    size_t step_count;
} TMConfig;

typedef enum {
    SC_L       = 0,
    SC_NL      = 1,
    SC_RL      = 2,
    SC_BPL     = 3,
    SC_PL      = 4,
    SC_PSPACE  = 5,
    SC_NPSPACE = 6,
    SC_EXPSPACE= 7,
    SC_SL      = 8,
    SC_CO_NL   = 9,
    SC_NC1     = 10,
    SC_NC2     = 11,
    SC_RNC     = 12,
    SC_BPP     = 13,
    SC_RP      = 14,
    SC_ZPP     = 15,
} SpaceClass;

typedef struct {
    size_t num_vertices;
    size_t num_edges;
    size_t *adj_offset;
    size_t *edge_targets;
    size_t start_vertex;
    size_t *accept_vertices;
    size_t num_accept_v;
} ConfigGraph;

typedef struct {
    size_t (*f)(size_t n);
    bool   (*verify)(size_t n);
    size_t space_usage;
} SpaceConstructible;

typedef struct {
    size_t n;
    size_t m;
    size_t k;
    size_t num_functions;
    bool  *hash_data;
} PairwiseHashFamily;

typedef struct {
    size_t seed_len;
    size_t output_len;
    size_t space_bound;
    double epsilon;
    bool  *generator_data;
    size_t data_size;
} SpacePRG;

typedef struct {
    size_t n;
    size_t d;
    size_t *edges;
    bool   directed;
    double *eigenvalues;
    bool   eigen_computed;
} RegularGraph;

typedef struct {
    double lambda1;
    double lambda2;
    double lambda_max;
    double spectral_gap;
    double edge_expansion;
    double mixing_rate;
} SpectralData;

typedef struct {
    size_t width;
    size_t layers;
    size_t *trans_0;
    size_t *trans_1;
    size_t start_node;
    size_t accept_node;
    size_t *accept_set;
    size_t num_accept;
    bool  *input_map;
    size_t input_count;
} BranchingProgram;

bool sd_savitch_simulate(const SpaceTM *tm, const char *input, size_t ilen, bool *result);
bool sd_immerman_szelepcsenyi_verify(const SpaceTM *tm, const char *input, size_t ilen, bool *result);
bool sd_reingold_ustconn(size_t n, const size_t *edges, size_t m, size_t s, size_t t, bool *connected);
const char *sd_class_name(SpaceClass sc);
bool sd_class_contains(SpaceClass outer, SpaceClass inner);
size_t sd_savitch_space_bound(size_t nspace_bound);
double sd_error_amplify(double base_error, size_t k_repetitions);
bool sd_verify_simulation(const SpaceTM *original, const SpaceTM *simulator, double epsilon);
ConfigGraph *sd_config_graph_build(const SpaceTM *tm, const char *input, size_t ilen);
void sd_config_graph_free(ConfigGraph *cg);
size_t sd_config_graph_reachable(const ConfigGraph *cg, size_t start, size_t steps, bool *reachable);
size_t sd_config_graph_bfs(const ConfigGraph *cg, size_t *out, size_t max_out);
bool sd_path_exists_dag(const size_t *adj_off, const size_t *adj_edges, size_t n, size_t s, size_t t, size_t depth, bool *result);
size_t sd_count_reachable_steps(const ConfigGraph *cg, size_t start, size_t t);
bool sd_verify_unreachable(const ConfigGraph *cg, size_t s, size_t w, size_t t, size_t reachable_count);
bool sd_nisan_prg_generate(const bool *seed, size_t slen, bool *out, size_t olen);
size_t sd_nisan_seed_len(size_t space, size_t output_len, double epsilon);
size_t sd_nisan_output_len(size_t seed_len, size_t space_bound);
bool sd_derandomize_rl(const SpaceTM *m, const char *input, size_t ilen, bool *result);
bool sd_derandomize_bpl(const SpaceTM *m, const char *input, size_t ilen, bool *result);
double sd_space_derand_overhead(size_t space, size_t time);
bool sd_log_constructible(size_t n);
size_t sd_ceil_log2(size_t n);
SpaceTM *sd_tm_create(size_t space, size_t num_states, size_t alpha_size);
void sd_tm_free(SpaceTM *tm);
bool sd_tm_set_transition(SpaceTM *tm, int state, char read_sym, int next_state, char write_sym, int direction);
void sd_config_print(const TMConfig *cfg);
bool sd_compute_eigenvalues(double *matrix, size_t n, double *eigenvalues, size_t max_iter, double tol);
PairwiseHashFamily *sd_hash_family_create(size_t n, size_t m);
size_t sd_hash_eval(const PairwiseHashFamily *phf, size_t idx, size_t x);
void sd_hash_family_free(PairwiseHashFamily *phf);
bool sd_verify_pairwise_independent(const PairwiseHashFamily *phf, size_t trials);
RegularGraph *sd_regular_graph_create(size_t n, size_t d, const size_t *edges);
void sd_regular_graph_free(RegularGraph *g);
SpectralData sd_compute_spectral(const RegularGraph *g);
bool sd_is_expander(const RegularGraph *g, double lambda_bound);
RegularGraph *sd_cycle_graph(size_t n);
RegularGraph *sd_complete_graph(size_t n);
RegularGraph *sd_zigzag_product(const RegularGraph *G, const RegularGraph *H);
RegularGraph *sd_graph_power(const RegularGraph *G, size_t t);
RegularGraph *sd_tensor_product(const RegularGraph *G, const RegularGraph *H);
RegularGraph *sd_replacement_product(const RegularGraph *G, const RegularGraph *H);
RegularGraph *sd_logspace_connectivity_network(size_t n);
bool sd_is_connected(const RegularGraph *G);
size_t sd_shortest_path(const RegularGraph *G, size_t s, size_t t);
size_t *sd_connected_components(const RegularGraph *G, size_t *num_comp);
bool sd_bp_random_walk(const BranchingProgram *bp, size_t start, size_t steps, const bool *rand_bits, size_t *final_node);
bool sd_bp_evaluate(const BranchingProgram *bp, const bool *input, size_t ilen);
BranchingProgram *sd_bp_create(size_t width, size_t layers);
void sd_bp_set_transition(BranchingProgram *bp, size_t layer, size_t from, bool bit, size_t to);
void sd_bp_free(BranchingProgram *bp);
double sd_bp_accept_probability(const BranchingProgram *bp, const double *input_probs, size_t ilen);
size_t sd_bp_diameter_bound(const BranchingProgram *bp);
bool sd_matrix_multiply(const double *A, const double *B, double *C, size_t n);
bool sd_power_iteration(const double *A, size_t n, double *eigenvector, double *eigenvalue, size_t max_iter, double tol);
double sd_rayleigh_quotient(const double *A, size_t n, const double *x);
double sd_vector_norm(const double *v, size_t n);
void sd_vector_normalize(double *v, size_t n);
void sd_vector_ones(double *v, size_t n);
void sd_vector_set(double *v, size_t n, double val);
double sd_vector_dot(const double *a, const double *b, size_t n);
void sd_graph_matvec(const RegularGraph *G, const double *x, double *y);
double sd_second_eigenvalue(const RegularGraph *G, size_t max_iter, double tol);
bool sd_graph_isomorphic(const RegularGraph *G, const RegularGraph *H);

#endif
