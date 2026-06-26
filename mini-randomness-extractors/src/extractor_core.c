
/**
 * src/extractor_core.c -- Core extractor definitions and operations
 *
 * L1: Extractor, WeakSource, SeededExt, TwoSourceExt typedef implementations
 * L2: Min-entropy H_inf, Shannon entropy H, Renyi entropy H_alpha,
 *     statistical distance Delta, Hellinger distance, KL divergence D_KL
 * L3: GF(2) field operations, probability simplex operations
 * L4: Leftover Hash Lemma bound, RTD lower bound, NZ seed bound,
 *     Chernoff inequality for min-entropy verification
 * L5: Generic extractor evaluation via Carter-Wegman universal hashing
 *
 * Reference theorems:
 *   LHL: Impagliazzo-Levin-Luby (1989) -- Delta <= 2^{-(k-m)/2}
 *   RTD: Radhakrishnan-Ta-Shma (2000) -- m <= k + 2*log(eps) + 2
 *   NZ:  Nisan-Zuckerman (1996) -- d >= log(n-k) + 2*log(1/eps) - O(1)
 *   Chernoff: Pr[|X_bar - mu| >= eps] <= 2*exp(-2*n*eps^2)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>

#include "extractor_core.h"

/*---------------------------------------------------------*/
#define MAX_N_BITS_EXACT  16
#define BITS_PER_WORD     (8 * sizeof(size_t))

/*---------------------------------------------------------*/
uint64_t ext_lcg_next(uint64_t *state) {
    *state = (*state * 6364136223846793005ULL + 1442695040888963407ULL);
    return *state;
}
double ext_lcg_double(uint64_t *state) {
    return (double)(ext_lcg_next(state) >> 11) / (1ULL << 53);
}

/*======== Section 1: Bit-vector utilities ==========*/

void ext_bits_pack(const bool *bits, size_t n, size_t *vec) {
    if (!bits || !vec || n == 0) return;
    size_t words = (n + BITS_PER_WORD - 1) / BITS_PER_WORD;
    memset(vec, 0, words * sizeof(size_t));
    for (size_t i = 0; i < n; i++) {
        if (bits[i])
            vec[i / BITS_PER_WORD] |= ((size_t)1) << (i % BITS_PER_WORD);
    }
}

void ext_bits_unpack(const size_t *vec, size_t n, bool *bits) {
    if (!vec || !bits || n == 0) return;
    for (size_t i = 0; i < n; i++)
        bits[i] = (vec[i / BITS_PER_WORD] >> (i % BITS_PER_WORD)) & 1;
}

size_t ext_popcount_bool(const bool *bits, size_t n) {
    if (!bits) return 0;
    size_t cnt = 0;
    for (size_t i = 0; i < n; i++) cnt += bits[i] ? 1 : 0;
    return cnt;
}

size_t ext_popcount(size_t x) {
    x = x - ((x >> 1) & (size_t)~(size_t)0 / 3);
    x = (x & (size_t)~(size_t)0 / 15 * 3)
        + ((x >> 2) & (size_t)~(size_t)0 / 15 * 3);
    x = (x + (x >> 4)) & (size_t)~(size_t)0 / 255 * 15;
    return (x * ((size_t)~(size_t)0 / 255)) >> (sizeof(size_t) - 1) * 8;
}

void ext_bits_copy(const bool *src, bool *dst, size_t n) {
    if (src && dst && n > 0) memcpy(dst, src, n * sizeof(bool));
}

void ext_bits_xor(const bool *a, const bool *b, bool *r, size_t n) {
    for (size_t i = 0; i < n; i++) r[i] = a[i] ^ b[i];
}

bool ext_bits_eq(const bool *a, const bool *b, size_t n) {
    for (size_t i = 0; i < n; i++) if (a[i] != b[i]) return false;
    return true;
}

void ext_bits_print(const bool *bits, size_t n, const char *label) {
    if (label) printf("%s: ", label);
    for (size_t i = n; i > 0; i--) printf("%c", bits[i-1] ? '1' : '0');
    printf(" (%zu bits)\n", n);
}

/*======== Section 2: Random bits and distribution sampling ==========*/

void ext_random_bits(bool *bits, size_t n, uint64_t *seed) {
    for (size_t i = 0; i < n; i++) bits[i] = (ext_lcg_double(seed) < 0.5);
}

void ext_sample_distribution(double bias, size_t n, size_t nsamples,
                              size_t *freq, uint64_t *seed) {
    if (!freq || n == 0 || nsamples == 0) return;
    if (n > MAX_N_BITS_EXACT) {
        fprintf(stderr, "ext_sample_distribution: n=%zu > max %d\n",
                n, MAX_N_BITS_EXACT);
        return;
    }
    size_t nvals = (size_t)1 << n;
    memset(freq, 0, nvals * sizeof(size_t));
    for (size_t s = 0; s < nsamples; s++) {
        size_t val = 0;
        for (size_t b = 0; b < n; b++)
            if (ext_lcg_double(seed) < bias) val |= ((size_t)1) << b;
        freq[val]++;
    }
}

void ext_sample_block_source(const double *biases, size_t n,
                              size_t nsamples, size_t *freq,
                              uint64_t *seed) {
    if (!biases || !freq || n == 0 || nsamples == 0) return;
    if (n > MAX_N_BITS_EXACT) {
        fprintf(stderr, "ext_sample_block_source: n=%zu > max %d\n",
                n, MAX_N_BITS_EXACT);
        return;
    }
    size_t nvals = (size_t)1 << n;
    memset(freq, 0, nvals * sizeof(size_t));
    for (size_t s = 0; s < nsamples; s++) {
        size_t val = 0;
        for (size_t b = 0; b < n; b++)
            if (ext_lcg_double(seed) < biases[b]) val |= ((size_t)1) << b;
        freq[val]++;
    }
}

/*======== Section 3: Min-entropy H_inf(X) = -log2(max_x Pr[X=x]) ==========*/

double ext_min_entropy_from_freq(const size_t *freq, size_t nvals,
                                  size_t nsamples) {
    if (!freq || nvals == 0 || nsamples == 0) return -1.0;
    size_t max_freq = 0;
    for (size_t i = 0; i < nvals; i++)
        if (freq[i] > max_freq) max_freq = freq[i];
    if (max_freq == 0) return -1.0;
    return log2((double)nsamples) - log2((double)max_freq);
}

double ext_min_entropy(const bool *dist, size_t nsamples, size_t n) {
    if (!dist || nsamples == 0 || n == 0) return -1.0;
    if (n > MAX_N_BITS_EXACT) return -1.0;
    size_t nvals = (size_t)1 << n;
    size_t *freq = (size_t*)calloc(nvals, sizeof(size_t));
    if (!freq) return -1.0;
    for (size_t s = 0; s < nsamples; s++) {
        size_t val = 0;
        for (size_t b = 0; b < n; b++)
            if (dist[s * n + b]) val |= ((size_t)1) << b;
        freq[val]++;
    }
    double result = ext_min_entropy_from_freq(freq, nvals, nsamples);
    free(freq);
    return result;
}

/*======== Section 4: Shannon entropy H(X) = -sum p_x log2(p_x) ==========*/

double ext_shannon_entropy(const size_t *freq, size_t nvals,
                            size_t nsamples) {
    if (!freq || nvals == 0 || nsamples == 0) return -1.0;
    double H = 0.0, inv_N = 1.0 / (double)nsamples;
    for (size_t i = 0; i < nvals; i++) {
        if (freq[i] > 0) {
            double p = (double)freq[i] * inv_N;
            H -= p * log2(p);
        }
    }
    return H;
}

/*======== Section 5: Renyi entropy H_alpha ==========*/

double ext_renyi_entropy(const size_t *freq, size_t nvals,
                          size_t nsamples, double alpha) {
    if (!freq || nvals == 0 || nsamples == 0) return -1.0;
    if (fabs(alpha - 1.0) < 1e-12)
        return ext_shannon_entropy(freq, nvals, nsamples);
    double sum = 0.0, inv_N = 1.0 / (double)nsamples;
    for (size_t i = 0; i < nvals; i++) {
        if (freq[i] > 0) {
            double p = (double)freq[i] * inv_N;
            sum += pow(p, alpha);
        }
    }
    if (sum <= 0) return -1.0;
    return log2(sum) / (1.0 - alpha);
}

/*======== Section 6: Statistical distance Delta(P,Q) = 0.5*|P-Q|_1 ==========*/

double ext_statistical_distance(const size_t *freq1, const size_t *freq2,
                                 size_t nvals, size_t total1,
                                 size_t total2) {
    if (!freq1 || !freq2 || nvals == 0 || total1 == 0 || total2 == 0)
        return -1.0;
    double delta = 0.0;
    double inv1 = 1.0 / (double)total1;
    double inv2 = 1.0 / (double)total2;
    for (size_t i = 0; i < nvals; i++) {
        double p1 = (double)freq1[i] * inv1;
        double p2 = (double)freq2[i] * inv2;
        delta += fabs(p1 - p2);
    }
    return 0.5 * delta;
}

/*======== Section 7: Hellinger distance ==========*/

double ext_hellinger_distance(const size_t *freq1, const size_t *freq2,
                               size_t nvals, size_t total1,
                               size_t total2) {
    if (!freq1 || !freq2 || nvals == 0 || total1 == 0 || total2 == 0)
        return -1.0;
    double sum_sq = 0.0;
    double inv1 = 1.0 / (double)total1;
    double inv2 = 1.0 / (double)total2;
    for (size_t i = 0; i < nvals; i++) {
        double p1 = (double)freq1[i] * inv1;
        double p2 = (double)freq2[i] * inv2;
        double diff = sqrt(p1) - sqrt(p2);
        sum_sq += diff * diff;
    }
    return sqrt(sum_sq) / sqrt(2.0);
}

/*======== Section 8: KL divergence D_KL(P||Q) ==========*/

double ext_kl_divergence(const size_t *freq1, const size_t *freq2,
                          size_t nvals, size_t total1, size_t total2) {
    if (!freq1 || !freq2 || nvals == 0 || total1 == 0 || total2 == 0)
        return -1.0;
    double dkl = 0.0;
    double inv1 = 1.0 / (double)total1;
    double inv2 = 1.0 / (double)total2;
    for (size_t i = 0; i < nvals; i++) {
        if (freq1[i] > 0) {
            double p = (double)freq1[i] * inv1;
            if (freq2[i] == 0) return INFINITY;
            double q = (double)freq2[i] * inv2;
            dkl += p * log2(p / q);
        }
    }
    return dkl;
}

/*======== Section 9: Extractor construction and evaluation ==========*/

Extractor ext_create(size_t n, size_t d, size_t m, double eps) {
    Extractor e = { n, d, m, eps };
    if (n == 0 || d > n || m == 0 || m > n || eps <= 0.0 || eps > 1.0)
        e.n = 0;
    return e;
}

bool ext_eval(const Extractor *e, const bool *source_bits,
              const bool *seed, bool *out) {
    if (!e || !source_bits || !seed || !out) return false;
    if (e->n == 0 || e->m == 0) return false;

    size_t x = 0, s = 0;
    for (size_t i = 0; i < e->n; i++)
        if (source_bits[i]) x |= ((size_t)1) << i;
    for (size_t i = 0; i < e->d && i < BITS_PER_WORD; i++)
        if (seed[i]) s |= ((size_t)1) << i;

    /* Smallest prime p > 2^n */
    size_t p = ((size_t)1 << e->n) + 1;
    if (p < 3) p = 3;
    while (1) {
        bool is_prime = true;
        size_t limit = (size_t)sqrt((double)p);
        for (size_t div = 2; div <= limit; div++)
            if (p % div == 0) { is_prime = false; break; }
        if (is_prime) break;
        p++;
    }

    size_t a = (s % (p - 1)) + 1;
    size_t b = (s / (p - 1)) % p;
    size_t hash_val = ((a * x + b) % p) % ((size_t)1 << e->m);

    for (size_t i = 0; i < e->m; i++)
        out[i] = (hash_val >> i) & 1;
    return true;
}

double ext_error_bound(const Extractor *e, size_t n, double k) {
    if (!e || n == 0 || k <= e->m) return 1.0;
    (void)n;
    double eps_lhl = pow(2.0, -(k - (double)e->m) / 2.0);
    double eps_rtd = pow(2.0, -((k - (double)e->m) - 2.0) / 2.0);
    return (eps_lhl < eps_rtd) ? eps_lhl : eps_rtd;
}

size_t ext_seed_len(size_t n, double min_entropy, double eps) {
    if (n == 0 || eps <= 0.0 || eps >= 1.0) return 0;
    double diff = (double)n - min_entropy;
    if (diff < 0) diff = 0;
    size_t log_nk = (size_t)ceil(log2(diff + 1.0));
    size_t log_eps = (size_t)ceil(2.0 * log2(1.0 / eps));
    size_t d = log_nk + log_eps + 3;
    return (d < 1) ? 1 : d;
}

size_t ext_output_len(size_t n, double min_entropy) {
    (void)n;
    if (min_entropy <= 0) return 0;
    double m_theo = min_entropy - min_entropy / 2.0 + 2.0;
    size_t m = (size_t)m_theo;
    return (m < 1) ? 1 : m;
}

bool ext_is_valid(const Extractor *e) {
    if (!e) return false;
    return (e->n > 0 && e->d <= e->n && e->m > 0 &&
            e->m <= e->n && e->eps > 0.0 && e->eps <= 1.0);
}

void ext_free(Extractor *e) {
    if (e) { e->n = e->d = e->m = 0; e->eps = 0.0; }
}

/*======== Section 10: Entropy deficit & verification ==========*/

double ext_entropy_deficit(size_t n, double min_entropy) {
    return (double)n - min_entropy;
}

bool ext_verify_min_entropy(const size_t *freq, size_t nvals,
                             size_t nsamples, double k, double delta) {
    if (!freq || nvals == 0 || nsamples == 0) return false;
    if (delta <= 0.0 || delta >= 1.0) return false;
    size_t max_freq = 0;
    for (size_t i = 0; i < nvals; i++)
        if (freq[i] > max_freq) max_freq = freq[i];
    double p_max_emp = (double)max_freq / (double)nsamples;
    double eps_chernoff = sqrt(log(2.0 / delta) / (2.0 * nsamples));
    double p_max_upper = p_max_emp + eps_chernoff;
    if (p_max_upper > 1.0) p_max_upper = 1.0;
    return (p_max_upper <= pow(2.0, -k));
}

/*======== Section 11: GF(2) field ops ==========*/

void ext_gf2_add(const bool *a, const bool *b, size_t n, bool *result) {
    for (size_t i = 0; i < n; i++) result[i] = a[i] ^ b[i];
}

/** GF(2^n) multiplication modulo x^n + x + 1 */
bool ext_gf2_mul(const bool *a, const bool *b, size_t n, bool *result) {
    if (!a || !b || !result || n == 0) return false;
    if (n > 32) return false;
    size_t max_deg = 2 * n;
    bool *prod = (bool*)calloc(max_deg, sizeof(bool));
    if (!prod) return false;
    for (size_t i = 0; i < n; i++) {
        if (!a[i]) continue;
        for (size_t j = 0; j < n; j++)
            if (b[j]) prod[i + j] ^= 1;
    }
    for (size_t i = max_deg - 1; i >= n; i--) {
        if (prod[i]) {
            prod[i] ^= 1;
            if (n > 1 && i >= n - 1) prod[i - n + 1] ^= 1;
            if (i == n) prod[0] ^= 1;
        }
    }
    memset(result, 0, n * sizeof(bool));
    for (size_t i = 0; i < n; i++) result[i] = prod[i];
    free(prod);
    return true;
}

/*======== Section 12: Conversion helpers ==========*/

size_t ext_bits_to_size(const bool *bits, size_t n) {
    size_t v = 0;
    for (size_t i = 0; i < n; i++) if (bits[i]) v |= ((size_t)1) << i;
    return v;
}

void ext_size_to_bits(size_t val, size_t n, bool *bits) {
    for (size_t i = 0; i < n; i++) bits[i] = (val >> i) & 1;
}

/*======== Section 13: Self-test ==========*/

int ext_core_self_test(void) {
    int f = 0;

    /* Bit pack roundtrip */
    { bool o[10]={1,0,1,1,0,0,1,0,1,0}; size_t v[2]={0}; bool r[10]={0};
      ext_bits_pack(o,10,v); ext_bits_unpack(v,10,r);
      for(int i=0;i<10;i++) if(o[i]!=r[i]) f++; }

    /* Popcount */
    if(ext_popcount(0xDEADBEEF)!=24) f++;

    /* Min-entropy uniform -> log2(N) */
    { size_t fr[8]={10,10,10,10,10,10,10,10};
      if(fabs(ext_min_entropy_from_freq(fr,8,80)-3.0)>0.01) f++; }

    /* Min-entropy deterministic -> 0 */
    { size_t fr[4]={100,0,0,0};
      if(fabs(ext_min_entropy_from_freq(fr,4,100))>0.01) f++; }

    /* Shannon uniform -> log2(4)=2 */
    { size_t fr[4]={25,25,25,25};
      if(fabs(ext_shannon_entropy(fr,4,100)-2.0)>0.01) f++; }

    /* Shannon deterministic -> 0 */
    { size_t fr[4]={100,0,0,0};
      if(fabs(ext_shannon_entropy(fr,4,100))>0.01) f++; }

    /* Renyi alpha=2 uniform -> 2 */
    { size_t fr[4]={25,25,25,25};
      if(fabs(ext_renyi_entropy(fr,4,100,2.0)-2.0)>0.01) f++; }

    /* Statistical distance identical -> 0 */
    { size_t a[4]={25,25,25,25},b[4]={25,25,25,25};
      if(fabs(ext_statistical_distance(a,b,4,100,100))>0.001) f++; }

    /* Statistical distance disjoint -> 1 */
    { size_t a[4]={100,0,0,0},b[4]={0,100,0,0};
      if(fabs(ext_statistical_distance(a,b,4,100,100)-1.0)>0.001) f++; }

    /* KL divergence same -> 0 */
    { size_t a[4]={25,25,25,25},b[4]={25,25,25,25};
      if(fabs(ext_kl_divergence(a,b,4,100,100))>0.001) f++; }

    /* Hellinger same -> 0 */
    { size_t a[4]={25,25,25,25},b[4]={25,25,25,25};
      if(fabs(ext_hellinger_distance(a,b,4,100,100))>0.001) f++; }

    /* Entropy deficit */
    { if(fabs(ext_entropy_deficit(8,3.0)-5.0)>0.001) f++; }

    /* Valid extractor */
    { Extractor e=ext_create(8,3,4,0.1);
      if(!ext_is_valid(&e)) f++; }

    /* Invalid: d > n */
    { Extractor bad=ext_create(8,10,4,0.1);
      if(ext_is_valid(&bad)) f++; }

    /* Invalid: eps <= 0 */
    { Extractor bad=ext_create(8,3,4,0.0);
      if(ext_is_valid(&bad)) f++; }

    /* ext_eval deterministic */
    { Extractor e=ext_create(4,4,2,0.1);
      bool src[4]={1,0,1,0},seed[4]={0,0,1,1},o1[2]={0},o2[2]={0};
      ext_eval(&e,src,seed,o1); ext_eval(&e,src,seed,o2);
      if(o1[0]!=o2[0]||o1[1]!=o2[1]) f++; }

    /* Seed length bound */
    { size_t d=ext_seed_len(8,4.0,0.1);
      if(d<5||d>20) f++; }

    /* Output length bound */
    { size_t m=ext_output_len(8,4.0);
      if(m<1||m>8) f++; }

    /* GF(2) add = XOR */
    { bool a[4]={1,0,1,1},b[4]={0,1,1,0},r[4];
      ext_gf2_add(a,b,4,r);
      if(r[0]!=1||r[1]!=1||r[2]!=0||r[3]!=1) f++; }

    /* bits<->size_t roundtrip */
    { bool bs[8]={0,1,0,1,0,0,1,1},bk[8]={0};
      size_t v=ext_bits_to_size(bs,8); ext_size_to_bits(v,8,bk);
      for(int i=0;i<8;i++) if(bs[i]!=bk[i]) f++; }

    /* Min-entropy verification: uniform passes
     * 80 uniform samples over 8 values -> H_inf=3.0.
     * Chernoff with 80 samples gives ~0.15 interval width,
     * so verifying k=1.5 (much less than 3.0) should pass. */
    { size_t fr[8]={10,10,10,10,10,10,10,10};
      if(!ext_verify_min_entropy(fr,8,80,1.5,0.05)) f++; }

    /* Min-entropy verification: biased fails
     * Max freq = 73/80 = 0.9125. Upper bound ~1.0 > 2^{-k} for any k. */
    { size_t fr[8]={73,1,1,1,1,1,1,1};
      if(ext_verify_min_entropy(fr,8,80,3.0,0.05)) f++; }

    /* Error bound monotonicity */
    { Extractor e=ext_create(4,4,2,0.1);
      double eh=ext_error_bound(&e,4,6.0);
      double el=ext_error_bound(&e,4,3.0);
      if(eh>el) f++; }

    return f;
}
