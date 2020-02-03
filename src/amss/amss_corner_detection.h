#ifndef amss_corner_detection_h__
#define amss_corner_detection_h__

#include "chain.h"

list_ptr amss_corner_detection (float **, long, long, float, long, long, float, float**);

list_ptr test_detection (char *, char *);
void draw_corners (float ** src, long nx, long ny, float ** mask);

#endif // amss_corner_detection_h__
