#ifndef least_squares_h__
#define least_squares_h__

#include "chain.h"

void error_percentile(list_ptr, float);
void solve_normal_equations (long, long, float **, float *, float *);

#endif /* least_squares_h__ */
