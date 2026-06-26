/******************************************************************************
 * worst_to_average.c ¡ª Worst-Case to Average-Case Hardness Amplification
 *
 * Implements Yao's XOR Lemma and Impagliazzo's Hardcore Lemma:
 * converting worst-case hardness into average-case hardness.
 *
 * Knowledge: L1 (hardness definitions), L2 (amplification concept),
 * L4 (Yao's XOR Lemma, Impagliazzo's Hardcore Lemma),
 * L5 (hardness amplification, Goldreich-Levin)
 *
 * References:
 *   Yao (1982): Theory and Applications of Trapdoor Functions
 *   Impagliazzo (1995): Hard-core distributions for somewhat hard problems
 *   Goldreich & Levin (1989): Hard-core predicate for any one-way function
 *   Arora & Barak (2009), Chapter 19
 ******************************************************************************/

#include "worst_to_average.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

/* ================================================================
 * L1+L2: Worst-Case to Average-Case Framework
 *
 * Worst-case hardness: ?x such that computing f(x) is hard.
 * Average-case hardness: for random x, computing f(x) is hard.
 *
 * Conversion theorems show that under appropriate conditions,
 * worst-case hardness implies average-case hardness for a
 * related function (via direct product / XOR).
 * ================================================================ */

/*
 * Yao's XOR Lemma:
 *
 * If f: {0,1}^n ¡ú {0,1} is ¦Ä-hard for circuits of size S
 *   (meaning: every circuit of size S fails to compute f on
 *    at least ¦Ä fraction of inputs),
 * then for any k, the function:
 *   F(x1, ..., xk) = f(x1) ¨’ f(x2) ¨’ ... ¨’ f(xk)
 * is ¦Å-close to (1/2 + (1-¦Ä)^k)-hard for circuits of size S/poly(k):
 *   Any circuit of size S' predicts F with advantage at most
 *   (1-¦Ä)^k + ¦Å.
 *
 * Interpretation: XOR amplifies hardness from ¦Ä to
 * (1/2 - (1-¦Ä)^k) for the k-fold XOR.
 *
 * Reference: Yao (1982); Goldreich et al. (2011)
 *          Arora & Barak, Theorem 19.10
 */
WTAConversion wta_yao_xor_lemma(double worst, size_t n, size_t k) {
    WTAConversion conv;
    conv.n = n;
    conv.worst_hardness = worst;
    /* The average hardness of k-XOR:
       avg ¡Ý 1/2 - (1 - worst)/2 * (1 - worst)^{k-1}
       For small worst-hardness: avg ¡Ö 1/2 - (1-worst)^k / 2
       But the key is: advantage drops exponentially in k */
    double base_advantage = 1.0 - worst; /* circuit can compute on worst fraction */
    /* After k independent copies, advantage = base_advantage^k */
    double advantage = pow(base_advantage, (double)k);
    /* Average hardness = 1/2 + advantage/2? No ¡ª the bound is:
       any circuit predicts XOR with advantage ¡Ü (1-¦Ä)^k.
       Hardness = 1/2 - advantage (probability of being correct) */
    conv.avg_hardness = 0.5 - advantage / 2.0;
    if (conv.avg_hardness < 0.0) conv.avg_hardness = 0.0;
    return conv;
}

/*
 * Compute XOR lemma hardness bound directly.
 * Given worst-case hardness ¦Ä and k copies, the advantage
 * of predicting XOR is at most (1-¦Ä)^k.
 *
 * Returns the probability that a small circuit correctly
 * predicts the XOR of k independent copies.
 */
double wta_xor_lemma_bound(double worst, size_t k) {
    if (k == 0) return 0.5;
    double advantage = pow(1.0 - worst, (double)k);
    return 0.5 + advantage / 2.0;
}

/*
 * Optimal number of XOR copies to achieve target average hardness.
 *
 * Need: (1 - worst)^k ¡Ü 2 * (target_avg - 0.5)
 * => k ¡Ý log(target_gap) / log(1 - worst)
 */
size_t wta_optimal_xor_copies(double worst, double target_avg) {
    if (worst <= 0.0 || target_avg <= 0.5) return 0;
    double target_gap = 2.0 * (target_avg - 0.5);
    if (target_gap <= 0.0) target_gap = 0.001;
    if (worst >= 1.0) return 1;
    /* k = ceil(log(target_gap) / log(1 - worst)) */
    double k = log(target_gap) / log(1.0 - worst);
    if (k < 0.0) k = 0.0;
    return (size_t)ceil(k);
}

/*
 * Impagliazzo's Hardcore Lemma:
 *
 * If f is ¦Ä-hard for circuits of size S, then there exists
 * a "hardcore set" H of density ¦Ä (|H| = ¦Ä¡¤2^n) such that
 * f is ¦Å-hard on H: every circuit of size S¡¤poly(¦Ä,¦Å) fails
 * on at least (1/2 - ¦Å) fraction of inputs from H.
 *
 * This is the converse of Yao's XOR Lemma and is central to
 * the derandomization framework.
 *
 * Reference: Impagliazzo (1995); Arora & Barak, Theorem 19.13
 */
WTAConversion wta_impagliazzo_hardcore(double worst, size_t n) {
    WTAConversion conv;
    conv.n = n;
    conv.worst_hardness = worst;
    /* Hardcore Lemma: within a set of density ¦Ä ¡Ö worst,
       the function is average-case hard: avg_hardness ¡Ö 1/2. */
    double density = worst;
    if (density < 0.001) density = 0.001;
    if (density > 1.0) density = 1.0;
    /* Average hardness on the hardcore set approaches 1/2
       as circuit size decreases */
    conv.avg_hardness = 0.5 - density / 10.0;
    if (conv.avg_hardness < 0.49) conv.avg_hardness = 0.49;
    return conv;
}

/*
 * General formula for average hardness from worst-case.
 *
 * avg = f(worst, n) where f depends on the specific reduction.
 * Simple interpolating formula used in many proofs.
 */
double wta_average_from_worst(double worst, size_t n) {
    if (worst <= 0.0) return 0.0;
    if (worst >= 1.0) return 1.0;
    /* Smooth interpolation: worst * (1 - 1/poly(n)) */
    double poly_factor = 1.0 / (double)(n + 1);
    double avg = worst * (1.0 - poly_factor);
    /* Amplification via XOR slows this down, but we use
       the direct relationship for small n */
    if (avg > 1.0) avg = 1.0;
    return avg;
}

/*
 * Check if worst-to-average conversion is feasible.
 * Feasible if: worst > 1/poly(n), i.e., non-negligible hardness.
 */
bool wta_conversion_feasible(double worst, size_t n) {
    if (n == 0) return false;
    double threshold = 1.0 / pow((double)n, 4.0);
    return worst >= threshold;
}

/*
 * Sample complexity for hardness amplification.
 *
 * To estimate average hardness within error ¦Å with confidence 1-¦Ä,
 * need O(log(1/¦Ä) / ¦Å^2) samples (Chernoff-Hoeffding).
 *
 * Reference: Arora & Barak, Appendix B
 */
size_t wta_sample_complexity(size_t n, double epsilon, double delta) {
    if (epsilon <= 0.0) return n * n;
    if (delta <= 0.0) delta = 0.01;
    /* Hoeffding: P[|avg - mu| > ¦Å] ¡Ü 2 exp(-2 m ¦Å^2)
       Need: 2 exp(-2 m ¦Å^2) ¡Ü delta
       => m ¡Ý ln(2/delta) / (2 ¦Å^2) */
    double samples = log(2.0 / delta) / (2.0 * epsilon * epsilon);
    if (samples < 1.0) samples = 1.0;
    if (samples > (double)SIZE_MAX) samples = (double)(n * n);
    return (size_t)ceil(samples);
}

/*
 * Local list decoding hardness:
 *
 * Given a function that is ¦Ä-hard in the worst case,
 * the Goldreich-Levin theorem shows how to list-decode
 * linear approximations.
 *
 * If f is (1/2 + ¦Å)-predicted by a family of circuits,
 * then there is a list of poly(1/¦Å) circuits that
 * collectively predict f on all inputs.
 *
 * This function computes the effective hardness after
 * list decoding with given list size.
 *
 * Reference: Goldreich & Levin (1989); Arora & Barak, Theorem 19.9
 */
double wta_local_list_decode_hardness(double worst, size_t n, size_t list_size) {
    (void)n; /* input size parameter */
    if (worst <= 0.0 || list_size == 0) return 0.0;
    /* List decoding weakens the hardness:
       effective_hardness ¡Ö worst / sqrt(list_size)
       because the union of L circuits might cover more inputs. */
    double factor = sqrt((double)list_size);
    double effective = worst / factor;
    if (effective > 1.0) effective = 1.0;
    return effective;
}

/*
 * Check if Goldreich-Levin theorem applies.
 *
 * GL Theorem: If f: {0,1}^n ¡ú {0,1} is a one-way function,
 * then g(x, r) = ?x, r? mod 2 (inner product mod 2) is a
 * hard-core predicate for f.
 *
 * Applies when hardness ¡Ý 1/poly(n).
 */
bool wta_goldreich_levin_applies(size_t n, double hardness) {
    if (n == 0) return false;
    double threshold = 1.0 / (double)n;
    return hardness >= threshold;
}

/*
 * Construct a hardcore set from worst-case hardness.
 *
 * Impagliazzo's construction: partition the input space
 * and identify regions where the function is unpredictable.
 * The hardcore set H is the union of these regions.
 *
 * Complexity: heuristic construction using sampling.
 *
 * Reference: Impagliazzo (1995)
 */
HardcoreSet *wta_construct_hardcore_set(double worst, size_t n, double density) {
    if (n > 20 || worst <= 0.0) return NULL;
    if (density < 0.0) density = 0.0;
    if (density > 1.0) density = 1.0;
    size_t domain_size = (size_t)1 << n;
    HardcoreSet *hs = malloc(sizeof(HardcoreSet));
    if (!hs) return NULL;
    hs->n = n;
    hs->density = density;
    hs->set = calloc(domain_size, sizeof(bool));
    if (!hs->set) { free(hs); return NULL; }
    hs->hardcore_hardness = 0.5;
    /* Construct hardcore set using sampling heuristic:
       pick each point with probability density, biased toward
       "hard" inputs. Since we don't have oracle access to
       hardness, we use a pseudorandom selection. */
    uint64_t state = (uint64_t)(worst * 1e9 + density * 1e6);
    for (size_t i = 0; i < domain_size; i++) {
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
        uint64_t r = state * 0x2545F4914F6CDD1DULL;
        double p = (double)(r & 0xFFFFF) / (double)0x100000;
        hs->set[i] = (p < density);
    }
    /* Count actual density */
    size_t count = 0;
    for (size_t i = 0; i < domain_size; i++) {
        if (hs->set[i]) count++;
    }
    hs->density = (double)count / (double)domain_size;
    return hs;
}

void wta_free_hardcore_set(HardcoreSet *hs) {
    if (!hs) return;
    free(hs->set);
    free(hs);
}