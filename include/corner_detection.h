#ifndef corner_detection_h__
#define corner_detection_h__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define pi 3.1415927


void dummy_mirror (float **v, long nx, long ny);

/* creates dummy boundaries by mirroring */

/* ---------------------------------------------------------------------- */

void gauss_conv (float sigma, long nx, long ny, float hx, float hy, float precision, long bc, float **f);

/*
 Gaussian convolution.
*/


/* ------------------------------------------------------------------------ */

void struct_tensor (
float ***v, long nc, long nx, long ny, float hx, float hy, float sigma, float rho, float **dxx, float **dxy, float **dyy);

/*
 Calculates the structure tensor.
*/

/*--------------------------------------------------------------------------*/

void corner_detection (
float ***f, long nc, long nx, long ny, float hx, float hy, long type, float T, float sigma, float rho, float **v);

/*
 calculates structure tensor based corner detector;
 v[i][j]=255 for corners;
 v[i][j]=0   else;
*/


int process_image (char *in, char *out, long type, float T, float sigma, float rho);

/*--------------------------------------------------------------------------*/

void draw_corners (float ***u, long nc, long nx, long ny, float **v);

#endif // corner_detection_h__