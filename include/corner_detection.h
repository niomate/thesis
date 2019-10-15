#ifndef corner_detection_h__
#define corner_detection_h__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define pi 3.1415927

void corner_detection(
        float ***f, long nc, long nx, long ny, float hx, float hy, long type, float T, float sigma, float rho,
        float **v);

void corner_detection_quantile(
        float ***f, long nc, long nx, long ny, float hx, float hy, long type, float quant, float sigma, float rho,
        float **v);

int detect_and_draw_corners(char *in, char *out, long type, float T, float sigma, float rho);

void draw_corners(float ***u, long nc, long nx, long ny, float **v);


#endif // corner_detection_h__