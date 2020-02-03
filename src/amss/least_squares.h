#ifndef least_squares_h__
#define least_squares_h__

#include "chain.h"
#include <stdlib.h>

void error_percentile(list_ptr, float);
void solve_normal_equations(long, long, float**, float*, float*);
int linear_regression(float* x, float* y, float* a, float* b, float* error, size_t n);

#endif /* least_squares_h__ */
