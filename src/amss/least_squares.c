#include "least_squares.h"
#include "../utils.h"
#include "chain.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define sqr(x) (x * x)

// http://mathworld.wolfram.com/LeastSquaresFitting.html
int linear_regression(float* tau, float* d, float* slope, float* intercept, float* r_sqr, size_t n)
{
    float sum_d = 0;  /* Sum of dn */
    float sum_t = 0;  /* Sum of tn */
    float sum_d2 = 0; /* Sum of dn^2 */
    float sum_t2 = 0; /* Sum of tn^2 */
    float sum_dt = 0; /* Sum of tn*dn */

    for (size_t i = 0; i < n; ++i) {
        sum_d += d[i];
        sum_t += tau[i];
        sum_d2 += sqr(d[i]);
        sum_t2 += sqr(tau[i]);
        sum_dt += d[i] * tau[i];
    }

    /* Mean values of d and tau */
    float mean_d = sum_d / n;
    float mean_t = sum_t / n;

    /* Compute sums of squares */
    float ssdd, sstt, ssdt;
    ssdd = sum_d2 - n * sqr(mean_d);
    sstt = sum_t2 - n * sqr(mean_t);
    ssdt = sum_dt - n * mean_d * mean_t;

    /* Singular matrix ? */
    if (sstt == 0) {
        return 0;
    }

    /* Compute the slope from sums of squares */
    *slope = ssdt / sstt;
    *intercept = mean_d - (*slope) * mean_t;

    *r_sqr = sqr(ssdt) / (ssdd * sstt);
    return 1;
}

/*--------------------------------------------------------------------------*/

int floatcomp_asc(const void* a, const void* b)
{
    if (*(const float*)a < *(const float*)b) {
        return -1;
    }
    return *(const float*)a > *(const float*)b;
}

int floatcomp_desc(const void* a, const void* b)
{
    return -floatcomp_asc(a, b);
}

void error_percentile(list_ptr chains, float percentile)
{
    assert(percentile >= 0);
    assert(percentile <= 1);

    if (chains->size <= 5 || percentile == 1) {
        return;
    }

    float error[chains->size];
    int i = 0;
    for (node_ptr current = list_head(chains); current != NULL; current = current->next) {
        error[i++] = current->error;
    }

    qsort(error, chains->size, sizeof(float), floatcomp_desc);

    long index = (long)(percentile * chains->size);
    float threshold = error[index];

    for (node_ptr current = list_head(chains); current != NULL; current = current->next) {
        if (current->error < threshold) {
            list_delete(chains, current);
        }
    }
}
