#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                    AFFINE MORPHOLOGICAL SCALE-SPACE                      */
/*                                                                          */
/*                  (Copyright Joachim Weickert, 1/2015)                    */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/*
 Explicit scheme with central spatial differences.
 Experimentally stable for time steps ht <= 0.1.
 No numerical guarantee for an extremum principle.
*/

/*--------------------------------------------------------------------------*/

void alloc_matrix

(float ***matrix, /* matrix */
 long n1,         /* size in direction 1 */
 long n2)         /* size in direction 2 */

/* allocates memory for matrix of size n1 * n2 */


{
    long i;

    *matrix = (float **)malloc (n1 * sizeof (float *));
    if (*matrix == NULL) {
        printf ("alloc_matrix: not enough memory available\n");
        exit (1);
    }
    for (i = 0; i < n1; i++) {
        (*matrix)[i] = (float *)malloc (n2 * sizeof (float));
        if ((*matrix)[i] == NULL) {
            printf ("alloc_matrix: not enough memory available\n");
            exit (1);
        }
    }
    return;
}

/*--------------------------------------------------------------------------*/

void disalloc_matrix

(float **matrix, /* matrix */
 long n1,        /* size in direction 1 */
 long n2)        /* size in direction 2 */

/* disallocates memory for matrix of size n1 * n2 */

{
    long i;

    for (i = 0; i < n1; i++)
        free (matrix[i]);

    free (matrix);

    return;
}

/*--------------------------------------------------------------------------*/

void read_string

(char *v) /* string to be read */

/*
 reads a long value v
*/

{
    fgets (v, 80, stdin);
    if (v[strlen (v) - 1] == '\n')
        v[strlen (v) - 1] = 0;
    return;
}

/*--------------------------------------------------------------------------*/

void read_long

(long *v) /* value to be read */

/*
 reads a long value v
*/

{
    char row[80]; /* string for reading data */

    fgets (row, 80, stdin);
    if (row[strlen (row) - 1] == '\n')
        row[strlen (row) - 1] = 0;
    sscanf (row, "%ld", &*v);
    return;
}

/*--------------------------------------------------------------------------*/

void read_float

(float *v) /* value to be read */

/*
 reads a float value v
*/

{
    char row[80]; /* string for reading data */

    fgets (row, 80, stdin);
    if (row[strlen (row) - 1] == '\n')
        row[strlen (row) - 1] = 0;
    sscanf (row, "%f", &*v);
    return;
}

/*--------------------------------------------------------------------------*/

void read_pgm_and_allocate_memory

(const char *file_name, /* name of pgm file */
 long *nx,              /* image size in x direction, output */
 long *ny,              /* image size in y direction, output */
 float ***u)            /* image, output */

/*
  reads a greyscale image that has been encoded in pgm format P5;
  allocates memory for the image u;
  adds boundary layers of size 1 such that
  - the relevant image pixels in x direction use the indices 1,...,nx
  - the relevant image pixels in y direction use the indices 1,...,ny
*/

{
    FILE *inimage; /* input file */
    char row[80];  /* for reading data */
    long i, j;     /* loop variables */

    /* open file */
    inimage = fopen (file_name, "rb");
    if (NULL == inimage) {
        printf ("could not open file '%s' for reading, aborting.\n", file_name);
        exit (1);
    }

    /* read header */
    fgets (row, 80, inimage); /* skip format definition */
    fgets (row, 80, inimage);
    while (row[0] == '#') /* skip comments */
        fgets (row, 80, inimage);
    sscanf (row, "%ld %ld", nx, ny); /* read image size */
    fgets (row, 80, inimage);        /* read maximum grey value */

    /* allocate memory */
    alloc_matrix (u, (*nx) + 2, (*ny) + 2);

    /* read image data row by row */
    for (j = 1; j <= (*ny); j++)
        for (i = 1; i <= (*nx); i++)
            (*u)[i][j] = (float)getc (inimage);

    /* close file */
    fclose (inimage);

    return;

} /* read_pgm_and_allocate_memory */

/*--------------------------------------------------------------------------*/

void comment_line

(char *comment,    /* comment string (output) */
 char *lineformat, /* format string for comment line */
 ...)              /* optional arguments */

/*
  Add a line to the comment string comment. The string line can contain plain
  text and format characters that are compatible with sprintf.
  Example call: print_comment_line(comment,"Text %f %d",float_var,int_var);
  If no line break is supplied at the end of the input string, it is added
  automatically.
*/

{
    char line[80];
    va_list arguments;

    /* get list of optional function arguments */
    va_start (arguments, lineformat);

    /* convert format string and arguments to plain text line string */
    vsprintf (line, lineformat, arguments);

    /* add line to total commentary string */
    strncat (comment, line, 80);

    /* add line break if input string does not end with one */
    if (line[strlen (line) - 1] != '\n')
        sprintf (comment, "%s\n", comment);

    /* close argument list */
    va_end (arguments);

    return;

} /* comment_line */

/*--------------------------------------------------------------------------*/

void write_pgm

(float **u,       /* image, unchanged */
 long nx,         /* image size in x direction */
 long ny,         /* image size in y direction */
 char *file_name, /* name of pgm file */
 char *comments)  /* comment string (set 0 for no comments) */

/*
  writes a greyscale image into a pgm P5 file;
*/

{
    FILE *outimage;     /* output file */
    long i, j;          /* loop variables */
    float aux;          /* auxiliary variable */
    unsigned char byte; /* for data conversion */

    /* open file */
    outimage = fopen (file_name, "wb");
    if (NULL == outimage) {
        printf ("Could not open file '%s' for writing, aborting\n", file_name);
        exit (1);
    }

    /* write header */
    fprintf (outimage, "P5\n"); /* format */
    if (comments != 0)
        fprintf (outimage, comments);        /* comments */
    fprintf (outimage, "%ld %ld\n", nx, ny); /* image size */
    fprintf (outimage, "255\n");             /* maximal value */

    /* write image data */
    for (j = 1; j <= ny; j++)
        for (i = 1; i <= nx; i++) {
            aux = u[i][j] + 0.499999; /* for correct rounding */
            if (aux < 0.0)
                byte = (unsigned char)(0.0);
            else if (aux > 255.0)
                byte = (unsigned char)(255.0);
            else
                byte = (unsigned char)(aux);
            fwrite (&byte, sizeof (unsigned char), 1, outimage);
        }

    /* close file */
    fclose (outimage);

    return;

} /* write_pgm */

/*--------------------------------------------------------------------------*/

void dummies

(float **u, /* image matrix */
 long nx,   /* size in x direction */
 long ny)   /* size in y direction */

/* creates dummy boundaries by mirroring */

{
    long i, j; /* loop variables */

    for (i = 1; i <= nx; i++) {
        u[i][0] = u[i][1];
        u[i][ny + 1] = u[i][ny];
    }

    for (j = 0; j <= ny + 1; j++) {
        u[0][j] = u[1][j];
        u[nx + 1][j] = u[nx][j];
    }

    return;
}

/*--------------------------------------------------------------------------*/

float sgn (float x)

/*
 sign function
*/

{
    float sign;

    if (x > 0.0)
        sign = 1.0;
    else if (x < 0)
        sign = -1.0;
    else
        sign = 0.0;

    return (sign);
}

/*--------------------------------------------------------------------------*/

void amss

(float ht,  /* time step size, 0 < ht <= 0.1 */
 long nx,   /* image dimension in x direction */
 long ny,   /* image dimension in y direction */
 float hx,  /* pixel size in x direction */
 float hy,  /* pixel size in y direction */
 float **u) /* input: original image ;  output: smoothed */

/*
 affine morphological scale-space, explicit scheme
*/

{
    long i, j;                   /* loop variables */
    float **f;                   /* u at old time level */
    float fx, fy, fxx, fxy, fyy; /* derivatives */
    float two_hx;                /* 2.0 * hx, time saver */
    float two_hy;                /* 2.0 * hx, time saver */
    float hx_sqr;                /* hx * hx, time saver */
    float hy_sqr;                /* hy * hy, time saver */
    float two_hx_hy;             /* 2.0 * hx * hy, time saver */
    float help;                  /* time saver */


    /* ---- allocate memory for f ---- */

    alloc_matrix (&f, nx + 2, ny + 2);


    /* ---- copy u into f ---- */

    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
            f[i][j] = u[i][j];


    /* ---- create reflecting dummy boundaries for f ---- */

    dummies (f, nx, ny);


    /* ---- loop ---- */

    /* calculate some time savers */
    two_hx = 2.0 * hx;
    two_hy = 2.0 * hy;
    hx_sqr = hx * hx;
    hy_sqr = hy * hy;
    two_hx_hy = 2.0 * hx * hy;

    /* loop */
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++) {
            /* central spatial derivatives */
            fx = (f[i + 1][j] - f[i - 1][j]) / two_hx;
            fy = (f[i][j + 1] - f[i][j - 1]) / two_hy;
            fxx = (f[i + 1][j] - 2.0 * f[i][j] + f[i - 1][j]) / hx_sqr;
            fyy = (f[i][j + 1] - 2.0 * f[i][j] + f[i][j - 1]) / hy_sqr;
            if (fx * fy < 0.0)
                fxy = (f[i + 1][j + 1] - f[i][j + 1] - f[i + 1][j] + f[i][j] + f[i - 1][j - 1] -
                       f[i][j - 1] - f[i - 1][j] + f[i][j]) /
                      two_hx_hy;
            else
                fxy = (-f[i - 1][j + 1] + f[i][j + 1] + f[i + 1][j] - f[i][j] - f[i + 1][j - 1] +
                       f[i][j - 1] + f[i - 1][j] - f[i][j]) /
                      two_hx_hy;

            /* evolution */
            help = fx * fx * fyy + fy * fy * fxx - 2.0 * fx * fy * fxy;
            u[i][j] = f[i][j] + ht * sgn (help) * pow (fabs (help), 0.33333333333);
        }


    /* ---- free memory for f ---- */

    disalloc_matrix (f, nx + 2, ny + 2);

    return;

} /* amss */

/*--------------------------------------------------------------------------*/

void analyse_grey

(float **u,   /* image, unchanged */
 long nx,     /* pixel number in x direction */
 long ny,     /* pixel number in y direction */
 float *min,  /* minimum, output */
 float *max,  /* maximum, output */
 float *mean, /* mean, output */
 float *std)  /* standard deviation, output */

/*
 computes minimum, maximum, mean, and standard deviation of a greyscale
 image u
*/

{
    long i, j;    /* loop variables */
    double help1; /* auxiliary variable */
    float help2;  /* auxiliary variable */

    *min = u[1][1];
    *max = u[1][1];
    help1 = 0.0;
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++) {
            if (u[i][j] < *min)
                *min = u[i][j];
            if (u[i][j] > *max)
                *max = u[i][j];
            help1 = help1 + (double)u[i][j];
        }
    *mean = (float)help1 / (nx * ny);

    *std = 0.0;
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++) {
            help2 = u[i][j] - *mean;
            *std = *std + help2 * help2;
        }
    *std = sqrt (*std / (nx * ny));

    return;

} /* analyse_grey */

/*--------------------------------------------------------------------------*/

int main ()

{
    char in[80];         /* for reading data */
    char out[80];        /* for reading data */
    float **u;           /* image */
    long k;              /* loop variable */
    long nx, ny;         /* image size in x, y direction */
    float ht;            /* time step size */
    long kmax;           /* largest iteration number */
    float max, min;      /* largest, smallest grey value */
    float mean;          /* average grey value */
    float std;           /* standard deviation */
    char comments[1600]; /* string for comments */


    printf ("\n");
    printf ("AFFINE MORPHOLOGICAL SCALE-SPACE, EXPLICIT SCHEME\n\n");
    printf ("**************************************************\n\n");
    printf ("    Copyright 2015 by Joachim Weickert            \n");
    printf ("    Dept. of Mathematics and Computer Science     \n");
    printf ("    Saarland University, Saarbruecken, Germany    \n\n");
    printf ("    All rights reserved. Unauthorised usage,      \n");
    printf ("    copying, hiring, and selling prohibited.      \n\n");
    printf ("    Send bug reports to                           \n");
    printf ("    weickert@mia.uni-saarland.de                  \n\n");
    printf ("**************************************************\n\n");


    /* ---- read input image (pgm format P5) ---- */

    printf ("input image (pgm):                ");
    read_string (in);
    read_pgm_and_allocate_memory (in, &nx, &ny, &u);


    /* ---- read parameters ---- */

    printf ("time step size (<=0.1):           ");
    read_float (&ht);

    printf ("number of iterations:             ");
    read_long (&kmax);

    printf ("output image (pgm):               ");
    read_string (out);
    printf ("\n");


    /* ---- analyse initial image ---- */

    analyse_grey (u, nx, ny, &min, &max, &mean, &std);
    printf ("initial image\n");
    printf ("minimum:          %8.2f \n", min);
    printf ("maximum:          %8.2f \n", max);
    printf ("mean:             %8.2f \n", mean);
    printf ("standard dev.:    %8.2f \n\n", std);


    /* ---- process image ---- */

    for (k = 1; k <= kmax; k++) {
        /* perform one iteration */
        amss (ht, nx, ny, 1.0, 1.0, u);

        /* analyse result */
        analyse_grey (u, nx, ny, &min, &max, &mean, &std);
        printf ("iteration number: %8ld \n", k);
        printf ("minimum:          %8.2f \n", min);
        printf ("maximum:          %8.2f \n", max);
        printf ("mean:             %8.2f \n", mean);
        printf ("standard dev.:    %8.2f \n\n", std);
    } /* for */


    /* ---- write output image (pgm format P5) ---- */

    /* write parameter values in comment string */
    comments[0] = '\0';
    comment_line (comments, "# AMSS, explicit scheme\n");
    comment_line (comments, "# initial image:  %s\n", in);
    comment_line (comments, "# ht:             %8.4f\n", ht);
    comment_line (comments, "# iterations:     %8ld\n", kmax);

    /* write image data */
    write_pgm (u, nx, ny, out, comments);
    printf ("output image %s successfully written\n\n", out);


    /* ---- free memory  ---- */

    disalloc_matrix (u, nx + 2, ny + 2);

    return (0);
}
