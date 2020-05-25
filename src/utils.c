#include "utils.h"
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define _USE_MATH_DEFINES

void alloc_vector(float **vector, long n) {
  *vector = (float *)malloc(n * sizeof(float));
  if (*vector == NULL) {
    printf("alloc_vector: not enough storage available\n");
    exit(1);
  }
  return;
}

void alloc_vector_long

    (long **vector, /* vector */
     long n)        /* size */

/* allocates storage for a vector of size n and type long */

{
  *vector = (long *)malloc(n * sizeof(long));
  if (*vector == NULL) {
    printf("alloc_vector_long: not enough storage available\n");
    exit(1);
  }
  return;
}

void alloc_matrix(float ***matrix, long nx, long ny) {
  long i;

  *matrix = (float **)malloc(nx * sizeof(float *));
  if (*matrix == NULL) {
    printf("alloc_matrix: not enough storage available\n");
    exit(1);
  }
  for (i = 0; i < nx; i++) {
    (*matrix)[i] = (float *)malloc(ny * sizeof(float));
    if ((*matrix)[i] == NULL) {
      printf("alloc_matrix: not enough storage available\n");
      exit(1);
    }
  }
  return;
}

void alloc_cubix(float ****cubix, long nx, long ny, long nz) {
  long i, j;

  *cubix = (float ***)malloc(nx * sizeof(float **));
  if (*cubix == NULL) {
    printf("alloc_cubix: not enough storage available\n");
    exit(1);
  }
  for (i = 0; i < nx; i++) {
    (*cubix)[i] = (float **)malloc(ny * sizeof(float *));
    if ((*cubix)[i] == NULL) {
      printf("alloc_cubix: not enough storage available\n");
      exit(1);
    }
    for (j = 0; j < ny; j++) {
      (*cubix)[i][j] = (float *)malloc(nz * sizeof(float));
      if ((*cubix)[i][j] == NULL) {
        printf("alloc_cubix: not enough storage available\n");
        exit(1);
      }
    }
  }
  return;
}

void disalloc_vector(float *vector, long n) {
  free(vector);
  return;
}

void disalloc_matrix(float **matrix, long nx, long ny) {
  long i;
  for (i = 0; i < nx; i++)
    free(matrix[i]);
  free(matrix);
  return;
}

void disalloc_cubix(float ***cubix, long nx, long ny, long nz) {
  long i, j;
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      free(cubix[i][j]);
  for (i = 0; i < nx; i++)
    free(cubix[i]);
  free(cubix);
  return;
}

/*--------------------------------------------------------------------------*/

void read_string

    (char *v) /* string to be read */

/*
 reads a long value v
*/

{
  fgets(v, 80, stdin);
  if (v[strlen(v) - 1] == '\n')
    v[strlen(v) - 1] = 0;
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

  fgets(row, 80, stdin);
  if (row[strlen(row) - 1] == '\n')
    row[strlen(row) - 1] = 0;
  sscanf(row, "%ld", &*v);
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

  fgets(row, 80, stdin);
  if (row[strlen(row) - 1] == '\n')
    row[strlen(row) - 1] = 0;
  sscanf(row, "%f", &*v);
  return;
}

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
  inimage = fopen(file_name, "rb");
  if (NULL == inimage) {
    printf("could not open file '%s' for reading, aborting.\n", file_name);
    exit(1);
  }

  /* read header */
  fgets(row, 80, inimage); /* skip format definition */
  fgets(row, 80, inimage);
  while (row[0] == '#') /* skip comments */
    fgets(row, 80, inimage);
  sscanf(row, "%ld %ld", nx, ny); /* read image size */
  fgets(row, 80, inimage);        /* read maximum grey value */

  /* allocate memory */
  alloc_matrix(u, (*nx) + 2, (*ny) + 2);

  /* read image data row by row */
  for (j = 1; j <= (*ny); j++)
    for (i = 1; i <= (*nx); i++)
      (*u)[i][j] = (float)getc(inimage);

  /* close file */
  fclose(inimage);

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
  va_start(arguments, lineformat);

  /* convert format string and arguments to plain text line string */
  vsprintf(line, lineformat, arguments);

  /* add line to total commentary string */
  strncat(comment, line, 80);

  /* add line break if input string does not end with one */
  if (line[strlen(line) - 1] != '\n')
    sprintf(comment, "%s\n", comment);

  /* close argument list */
  va_end(arguments);

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
  outimage = fopen(file_name, "wb");
  if (NULL == outimage) {
    printf("Could not open file '%s' for writing, aborting\n", file_name);
    exit(1);
  }

  /* write header */
  fprintf(outimage, "P5\n"); /* format */
  if (comments != 0)
    fprintf(outimage, comments);          /* comments */
  fprintf(outimage, "%ld %ld\n", nx, ny); /* image size */
  fprintf(outimage, "255\n");             /* maximal value */

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
      fwrite(&byte, sizeof(unsigned char), 1, outimage);
    }

  /* close file */
  fclose(outimage);

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

float sgn(float x)

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
  *std = sqrt(*std / (nx * ny));

  return;

} /* analyse_grey */

long clamp(long n, long lo, long hi) {
  if (n > hi)
    return hi;
  if (n < lo)
    return lo;
  return n;
}

void imgcpy(float **src, float **dest, size_t nx, size_t ny) {
  for (size_t i = 0; i < nx; ++i) {
    for (size_t j = 0; j < ny; ++j) {
      dest[i][j] = src[i][j];
    }
  }
}

int in_circle(long x, long y, long cx, long cy, float radius) {
  return powf(x - cx, 2.0f) + powf(y - cy, 2.0f) <= radius * radius;
}

int out_of_bounds(long i, long nx) {
  return i < 1 || i > nx;
}
