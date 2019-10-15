#include "utils.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


void alloc_vector (float **vector, long n) {
    *vector = (float *)malloc (n * sizeof (float));
    if (*vector == NULL) {
        printf ("alloc_vector: not enough storage available\n");
        exit (1);
    }
    return;
}


void alloc_matrix (float ***matrix, long nx, long ny) {
    long i;

    *matrix = (float **)malloc (nx * sizeof (float *));
    if (*matrix == NULL) {
        printf ("alloc_matrix: not enough storage available\n");
        exit (1);
    }
    for (i = 0; i < nx; i++) {
        (*matrix)[i] = (float *)malloc (ny * sizeof (float));
        if ((*matrix)[i] == NULL) {
            printf ("alloc_matrix: not enough storage available\n");
            exit (1);
        }
    }
    return;
}

void alloc_cubix (float ****cubix, long nx, long ny, long nz) {
    long i, j;

    *cubix = (float ***)malloc (nx * sizeof (float **));
    if (*cubix == NULL) {
        printf ("alloc_cubix: not enough storage available\n");
        exit (1);
    }
    for (i = 0; i < nx; i++) {
        (*cubix)[i] = (float **)malloc (ny * sizeof (float *));
        if ((*cubix)[i] == NULL) {
            printf ("alloc_cubix: not enough storage available\n");
            exit (1);
        }
        for (j = 0; j < ny; j++) {
            (*cubix)[i][j] = (float *)malloc (nz * sizeof (float));
            if ((*cubix)[i][j] == NULL) {
                printf ("alloc_cubix: not enough storage available\n");
                exit (1);
            }
        }
    }
    return;
}

void disalloc_vector (float *vector, long n) {
    free (vector);
    return;
}

void disalloc_matrix (float **matrix, long nx, long ny) {
    long i;
    for (i = 0; i < nx; i++)
        free (matrix[i]);
    free (matrix);
    return;
}

void disalloc_cubix (float ***cubix, long nx, long ny, long nz) {
    long i, j;
    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
            free (cubix[i][j]);
    for (i = 0; i < nx; i++)
        free (cubix[i]);
    free (cubix);
    return;
}

int cmpfunc (const void *a, const void *b) {
    return (*(float *)a - *(float *)b);
}

float quantile (float *arr, float q, long n) {
    qsort (arr, n, sizeof (float), cmpfunc);
    long qindex = (long)(q * n);
    printf ("Quantile index is %ld.\n", qindex);
    return arr[qindex];
}

/* Calculates a threshold s.t. <q> percent of all pixels are below this threshold
a.k.a. keep the top 1-q percent corners */
float quantile_2D (float **u, float q, long nx, long ny) {
    float *pixels = malloc (nx * ny * sizeof (float));
    /* Flatten image for quantile computation */
    for (long i = 0; i < nx; ++i) {
        for (long j = 0; j < ny; ++j) {
            pixels[i * ny + j] = u[i][j];
        }
    }
    return quantile (pixels, q, nx * ny);
}


void read_pgm_and_allocate_memory (char *filename, long *nx, long *ny, float ***u) {
    /* open pgm file and read header */
    char buf[80];

    FILE *image = fopen (filename, "r");

    if (image == NULL) {
        fprintf (stderr, "Image %s not found!", filename);
    }

    fgets (buf, 300, image);
    fgets (buf, 300, image);
    while (buf[0] == '#')
        fgets (buf, 300, image);
    sscanf (buf, "%ld %ld", nx, ny);
    fgets (buf, 300, image);

    /* allocate storage */
    alloc_matrix (u, *nx, *ny);

    /* read image data */
    for (long i = 0; i < *nx; i++)
        for (long j = 0; j < *ny; j++)
            (*u)[i][j] = (float)getc (image);

    fclose (image);
}


float MSE (float **u, float **v, long nx, long ny) {
    float acc = 0;
    for (long i = 0; i < nx; ++i) {
        for (long j = 0; j < ny; ++j) {
            acc += powf (u[i][j] - v[i][j], 2.0);
        }
    }
    return acc / (nx * ny);
}

float MSE_filenames (char *im1, char *im2) {
    float **u;
    float **v;
    long unx, uny, vnx, vny;

    read_pgm_and_allocate_memory (im1, &unx, &uny, &u);
    read_pgm_and_allocate_memory (im2, &vnx, &vny, &v);

    assert (unx == vnx && uny == vny);

    float mse = MSE (u, v, unx, uny);

    disalloc_matrix (u, unx, uny);
    disalloc_matrix (v, vnx, vny);

    return mse;
}