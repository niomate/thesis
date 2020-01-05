#include "../mask/main.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main ()

{
    char row[80];             /* for reading data */
    char in[80];              /* for reading data */
    char out[80];             /* for reading data */
    float ***u;               /* input image */
    float **v;                /* image with corner location */
    long i, j, m;             /* loop variables */
    long nx, ny;              /* image size in x, y direction */
    long nc;                  /* number of channels */
    long type;                /* type of corner detector */
    long count;               /* number of corners */
    float T;                  /* threshold */
    float sigma;              /* noise scale */
    float rho;                /* integration scale */
    unsigned char byte;       /* for data conversion */
    FILE *inimage, *outimage; /* input file, output file */

    printf ("\n");
    printf ("CORNER DETECTION WITH THE STRUCTURE TENSOR\n\n");

    /* ---- read input image (pgm format P5 or ppm format P6) ---- */

    /* read image name */
    printf ("input image (pgm, ppm):           ");
    gets (in);

    /* open pgm file and read header */
    inimage = fopen (in, "r");
    fgets (row, 300, inimage); /* image type: P5 or P6 */
    if ((row[0] == 'P') && (row[1] == '5'))
        nc = 1; /* P5: grey scale image */
    else if ((row[0] == 'P') && (row[1] == '6'))
        nc = 3; /* P6: colour image */
    else {
        printf ("unknown image format");
        exit (0);
    }
    fgets (row, 300, inimage); /* first comment row */
    while (row[0] == '#')
        fgets (row, 300, inimage);     /* further comment rows */
    sscanf (row, "%ld %ld", &nx, &ny); /* image size */
    fgets (row, 300, inimage);         /* # grey values */

    /* allocate storage */
    alloc_cubix (&u, nc, nx + 2, ny + 2);

    /* read image data */
    for (j = 1; j <= ny; j++)
        for (i = 1; i <= nx; i++)
            for (m = 0; m <= nc - 1; m++)
                u[m][i][j] = (float)getc (inimage);
    fclose (inimage);

    /* ---- read other parameters ---- */

    printf ("corner detector:\n");
    printf ("  (0) Rohr:             det(J)\n");
    printf ("  (1) Tomasi-Kanade:    lambda_2\n");
    printf ("  (2) Foerstner-Harris: det(J)/tr(J)\n");
    printf ("your choice:                      ");
    gets (row);
    sscanf (row, "%ld", &type);
    printf ("threshold T (>=0):                ");
    gets (row);
    sscanf (row, "%f", &T);
    printf ("noise scale sigma (>=0):          ");
    gets (row);
    sscanf (row, "%f", &sigma);
    printf ("integration scale rho (>=0):      ");
    gets (row);
    sscanf (row, "%f", &rho);
    if (nc == 1)
        printf ("output image (pgm):               ");
    else if (nc == 3)
        printf ("output image (ppm):               ");
    gets (out);
    printf ("\n");

    /* ---- process image ---- */

    /* allocate storage for corner detector */
    alloc_matrix (&v, nx + 2, ny + 2);

    /* corner detection */
    corner_detection (u, nc, nx, ny, 1.0, 1.0, type, T, sigma, rho, v);

    /* insert corner marks into the image */
    draw_corners (u, nc, nx, ny, v);

    /* count corners */
    count = 0;
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
            if (v[i][j] == 255.0)
                count = count + 1;
    printf ("number of corners:   %3ld\n", count);

    /* ---- write output image ---- */

    /* open file and write header (incl. filter parameters) */
    outimage = fopen (out, "w");
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
    fprintf (outimage, "# initial image:  %s\n", in);
    fprintf (outimage, "# sigma:          %8.4f\n", sigma);
    fprintf (outimage, "# rho:            %8.4f\n", rho);
    fprintf (outimage, "# T:              %8.4f\n", T);
    fprintf (outimage, "# corners:        %8ld\n", count);
    fprintf (outimage, "%ld %ld \n255\n", nx, ny);

    /* write image data and close file */
    for (j = 1; j <= ny; j++)
        for (i = 1; i <= nx; i++)
            for (m = 0; m <= nc - 1; m++) {
                if (u[m][i][j] < 0.0)
                    byte = (unsigned char)(0.0);
                else if (u[m][i][j] > 255.0)
                    byte = (unsigned char)(255.0);
                else
                    byte = (unsigned char)(u[m][i][j]);
                fwrite (&byte, sizeof (unsigned char), 1, outimage);
            }
    fclose (outimage);
    printf ("output image %s successfully written\n\n", out);

    /* ---- disallocate storage ---- */

    disalloc_cubix (u, nc, nx + 2, ny + 2);
    disalloc_matrix (v, nx + 2, ny + 2);

    return (0);
}