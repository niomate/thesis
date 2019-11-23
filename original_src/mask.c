#include "corner_detection.h"
#include "mask.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

float calc_compression_mask (char *in, char *out, long type, float T, float sigma, float rho, long size) {
    /*
    First: read input image
    Consider what parameters to use for corner detection
    Calculate corner map.
    */
    char buf[80];
    int nc;
    long nx, ny;
    long i, j, m;
    long count;
    float percentage;
    unsigned char byte;
    float ***u;
    float **v;
    float ***mask;

    printf ("CALCULATION OF COMPRESSION MASK USING STRUCTURE TENSOR BASED CORNER DETECTION\n\n");

    FILE *inimage = fopen (in, "r");
    fgets (buf, 300, inimage); /* image type: P5 or P6 */
    if ((buf[0] == 'P') && (buf[1] == '5'))
        nc = 1; /* P5: grey scale image */
    else if ((buf[0] == 'P') && (buf[1] == '6'))
        nc = 3; /* P6: colour image */
    else {
        printf ("unknown image format");
        exit (0);
    }
    fgets (buf, 300, inimage); /* first comment buf */
    while (buf[0] == '#')
        fgets (buf, 300, inimage);     /* further comment bufs */
    sscanf (buf, "%ld %ld", &nx, &ny); /* image size */
    fgets (buf, 300, inimage);         /* # grey values */

    /* allocate storage */
    alloc_cubix (&u, nc, nx + 2, ny + 2);

    /* read image data */
    for (j = 1; j <= ny; j++)
        for (i = 1; i <= nx; i++)
            for (m = 0; m <= nc - 1; m++)
                u[m][i][j] = (float)getc (inimage);
    fclose (inimage);

    alloc_matrix (&v, nx + 2, ny + 2);

    /* corner detection */
    corner_detection (u, nc, nx, ny, 1.0, 1.0, type, T, sigma, rho, v);

    count = 0;
    for (i = 1; i <= nx; ++i) {
        for (j = 1; j <= ny; ++j) {
            if (v[i][j] == 255.0) {
                ++count;
            }
        }
    }

    printf ("Found %ld corners and junctions.\n\n", count);

    /* Keep areas around detected corners */

    /* Initialise mask */
    alloc_cubix (&mask, nc, nx + 2, ny + 2);

    for (i = 0; i <= nx + 1; ++i) {
        for (j = 0; j <= ny + 1; ++j) {
            for (m = 0; m <= nc - 1; ++m) {
                mask[m][i][j] = -1.0;
            }
        }
    }

    long k, l;

    for (m = 0; m <= nc - 1; m++) {
        for (i = 5; i <= nx - 4; i++) {
            for (j = 5; j <= ny - 4; j++) {
                /* keep area around points of interest */
                if (v[i][j] == 255.0) {
                    if (size == 0) {
                        mask[m][i][j] = 255.0;
                    } else
                        for (k = i - 4; k <= i + 4; k++) {
                            for (l = j - 4; l <= j + 4; l++) {
                                if ((k - i) * (k - i) + (l - j) * (l - j) <= size)
                                    mask[m][k][l] = 255.0;
                            }
                        }
                }
            }
        }
    }

    /* calculate percentage of pixels kept */
    count = 0;
    for (i = 1; i <= nx; ++i) {
        for (j = 1; j <= ny; ++j) {
            if (nc == 3) {
                if (mask[0][i][j] != -1.0 || mask[1][i][j] != -1.0 || mask[2][i][j] != -1.0) {
                    ++count;
                }
            } else if (nc == 1) {
                if (mask[0][i][j] != -1.0) {
                    ++count;
                }
            }
        }
    }

    percentage = (100 * (float)count) / (nx * ny);

    printf ("Percentage of pixels kept: %3.2f\n", percentage);

    /* open file and write header (incl. filter parameters) */
    FILE *outimage = fopen (out, "w");
    if (nc == 1)
        fprintf (outimage, "P5 \n");
    else if (nc == 3)
        fprintf (outimage, "P6 \n");
    fprintf (outimage, "# structure tensor analysis\n");
    if (type == 0)
        fprintf (outimage, "# Rohr corner detector \n");
    if (type == 1)
        fprintf (outimage, "# Tomasi-Kanade corner detector \n");
    if (type == 2)
        fprintf (outimage, "# Foerstner-Harris corner detector \n");
    fprintf (outimage, "# initial image:     %s\n", in);
    fprintf (outimage, "# sigma:             %8.4f\n", sigma);
    fprintf (outimage, "# rho:               %8.4f\n", rho);
    fprintf (outimage, "# T:                 %8.4f\n", T);
    fprintf (outimage, "# corners:           %8ld\n", count);
    fprintf (outimage, "# percentage kept:   %8.4f\n", percentage);
    fprintf (outimage, "%ld %ld \n255\n", nx, ny);

    /* write image data and close file */
    for (j = 1; j <= ny; j++)
        for (i = 1; i <= nx; i++)
            for (m = 0; m <= nc - 1; m++) {
                if (mask[m][i][j] < 0.0)
                    byte = (unsigned char)(0.0);
                else if (mask[m][i][j] > 255.0)
                    byte = (unsigned char)(255.0);
                else
                    byte = (unsigned char)(mask[m][i][j]);
                fwrite (&byte, sizeof (unsigned char), 1, outimage);
            }
    fclose (outimage);
    printf ("output image %s successfully written\n\n", out);

    /* ---- disallocate storage ---- */

    disalloc_cubix (u, nc, nx + 2, ny + 2);
    disalloc_cubix (mask, nc, nx + 2, ny + 2);
    disalloc_matrix (v, nx + 2, ny + 2);

    return percentage;
}


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