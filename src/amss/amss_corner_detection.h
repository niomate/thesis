#ifndef amss_corner_detection_h__
#define amss_corner_detection_h__

#include "chain.h"

list_ptr amss_corner_detection (float **u, long nx, long ny, float ht, long t_0, long t_max, long step, float q, float **v);

list_ptr test_detection (char *in, char *out);

#endif // amss_corner_detection_h__
