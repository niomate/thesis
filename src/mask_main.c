#include "corner_detection.h"
#include "mask.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main () {
    char buf[80];
    char in[80];
    char out[80];
    long type;
    float rho, sigma, T;

    printf ("CORNER AND JUNCTION BASED COMPRESSION WITH EED INPAINTING");

    /* Read input image */
    printf ("\n\nInput image (pgm):\t");
    gets (in);

    printf ("corner detector:\n");
    printf ("  (0) Rohr:             det(J)\n");
    printf ("  (1) Tomasi-Kanade:    lambda_2\n");
    printf ("  (2) Foerstner-Harris: det(J)/tr(J)\n");
    printf ("your choice:                      ");
    gets (buf);
    sscanf (buf, "%ld", &type);
    printf ("threshold T (>=0):                ");
    gets (buf);
    sscanf (buf, "%f", &T);
    printf ("noise scale sigma (>=0):          ");
    gets (buf);
    sscanf (buf, "%f", &sigma);
    printf ("integration scale rho (>=0):      ");
    gets (buf);
    sscanf (buf, "%f", &rho);
    printf ("output image (pgm):               ");
    gets (out);
    printf ("\n");

    calc_compression_mask (in, out, type, T, sigma, rho, 5);

    return 0;
}