#include "utils.h"
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define _USE_MATH_DEFINES

/**
 * \brief Allocate a float array of size n
 *
 * @param vector: pointer to the float array
 * @param n: size of the array to be located
 */
void alloc_vector(float **vector, long n) {
  *vector = (float *)malloc(n * sizeof(float));
  if (*vector == NULL) {
    printf("alloc_vector: not enough storage available\n");
    exit(1);
  }
  return;
}

/**
 * \brief Allocate a long array of size n
 *
 * @param vector: pointer to the long array
 * @param n: size of the array to be located
 */
void alloc_vector_long(long **vector, long n) {
  *vector = (long *)malloc(n * sizeof(long));
  if (*vector == NULL) {
    printf("alloc_vector_long: not enough storage available\n");
    exit(1);
  }
  return;
}

/**
 * \brief Allocate a float matrix of size nx*ny
 *
 * @param vector: pointer to the matrix
 * @param nx: size in x dimension
 * @param ny: size in y dimension
 */
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

/**
 * \brief Allocate a float cubix of size nx*ny*nz
 *
 * @param vector: pointer to the cubix
 * @param nx: size in x dimension
 * @param ny: size in y dimension
 * @param nz: size in z dimension
 */
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

/**
 * \brief Free an array of size n starting at position of pointer vector
 *
 * @param vector: start of the array
 * @param n: size of the array
 */
void disalloc_vector(float *vector, long n) {
  free(vector);
  return;
}

/**
 * \brief Free a matrix of size n starting at position of pointer matrix
 *
 * @param vector: start of the array
 * @param nx: size in x dimension
 * @param ny: size in y dimension
 */
void disalloc_matrix(float **matrix, long nx, long ny) {
  long i;
  for (i = 0; i < nx; i++)
    free(matrix[i]);
  free(matrix);
  return;
}

/**
 * \brief Free a matrix of size n starting at position of pointer cubix
 *
 * @param vector: start of the array
 * @param nx: size in x dimension
 * @param ny: size in y dimension
 * @param nz: size in z dimension
 */
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

/**
 * \brief Read a string of size 80 into v from stdin
 *
 * @param v: buffer that is read into from stdin
 */
void read_string(char *v) {
  fgets(v, 80, stdin);
  if (v[strlen(v) - 1] == '\n')
    v[strlen(v) - 1] = 0;
}

/**
 * \brief Read long from stdin
 *
 * @param v: afterwards, contains number from stdin
 */
void read_long(long *v) {
  char row[80];
  fgets(row, 80, stdin);
  if (row[strlen(row) - 1] == '\n')
    row[strlen(row) - 1] = 0;
  sscanf(row, "%ld", &*v);
}

/**
 * \brief Read float from stdin
 *
 * @param v: afterwards, contains number from stdin
 */
void read_float(float *v) {
  char row[80];
  fgets(row, 80, stdin);
  if (row[strlen(row) - 1] == '\n')
    row[strlen(row) - 1] = 0;
  sscanf(row, "%f", &*v);
}

/**
 * \brief Read a pgm image and store it in u
 *
 * Reads a greyscale image that has been encoded in pgm format P5 and
 * allocates memory for the image u.
 * Afterwards, adds boundary layers of size 1 such that
 * the relevant image pixels in x direction use the indices 1,...,nx
 * and the relevant image pixels in y direction use the indices 1,...,ny.
 *
 * @param file_name: name of the image file
 * @param nx: size of the image in x dimension
 * @param ny: size of the image in y dimension
 * @param u: pointer to the image buffer
 */
void read_pgm_and_allocate_memory(const char *file_name, long *nx, long *ny,
                                  float ***u) {
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
}

/**
 * \brief Add a line to the comment string.
 *
 * Add a line to the comment string comment. The string line can contain plain
 * text and format characters that are compatible with sprintf.
 * Example call: print_comment_line(comment,"Text %f %d",float_var,int_var);
 * If no line break is supplied at the end of the input string, it is added
 * automatically.
 *
 * @param comment: comment string
 * @param lineform: format string
 * @param ...: optional arguments
 */
void comment_line(char *comment, char *lineformat, ...) {
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
}

/**
 * \brief Write image buffer to file.
 *
 * Write the given image buffer to the given file and include comments
 *
 * @param u: image buffer, unchanged
 * @param nx: size of the image in x dimension
 * @param ny: size of the image in y dimension
 * @param file_name: name of the image file
 * @param comments: comment string, NULL for no comments
 */
void write_pgm(float **u, long nx, long ny, char *file_name, char *comments) {
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
  if (comments != 0) {
    fprintf(outimage, comments); /* comments */
  }
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
}

/**
 * \brief Create reflecting boundary conditions.
 *
 * @param u: image buffer
 * @param nx: size of the image in x dimension
 * @param ny: size of the image in y dimension
 */
void dummies(float **u, long nx, long ny) {
  long i, j; /* loop variables */

  for (i = 1; i <= nx; i++) {
    u[i][0] = u[i][1];
    u[i][ny + 1] = u[i][ny];
  }

  for (j = 0; j <= ny + 1; j++) {
    u[0][j] = u[1][j];
    u[nx + 1][j] = u[nx][j];
  }
}

/**
 * \brief Sign function. -1 if negative, 1 if positive, 0 if 0
 *
 * @param x: number which sign should be calculated
 */
float sgn(float x) {
  float sign;
  if (x > 0.0)
    sign = 1.0;
  else if (x < 0)
    sign = -1.0;
  else
    sign = 0.0;
  return sign;
}

/**
 * \brief Compute statistics for the given image.
 *
 *  Compute minimum, maximum, mean and standard deviation of the greyscale
 *  image u.
 *
 * @param u: image buffer
 * @param nx: size of the image in x dimension
 * @param ny: size of the image in y dimension
 * @param min: minimum of u, output
 * @param max: maximum of u, output
 * @param mean: mean of u, output
 * @param std: standard deviation, output
 */
void analyse_grey(float **u, long nx, long ny, float *min, float *max,
                  float *mean, float *std) {
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
}

/**
 * \brief Clamp the given number between lo and hi
 *
 * Returns lo if n < lo, hi if n > hi and n else
 *
 * @param n: number to be clamped
 * @param lo: lower bound
 * @param hi: higher bound
 */
long clamp(long n, long lo, long hi) {
  if (n > hi)
    return hi;
  if (n < lo)
    return lo;
  return n;
}

/**
 * \brief Copy image from src to dest
 *
 * @param src: source image
 * @param dest: destination buffer
 * @param nx: size of the image in x dimension
 * @param ny: size of the image in y dimension
 */
void imgcpy(float **src, float **dest, size_t nx, size_t ny) {
  for (size_t i = 0; i < nx; ++i) {
    for (size_t j = 0; j < ny; ++j) {
      dest[i][j] = src[i][j];
    }
  }
}

/**
 * \brief Check if (x,y) is in a circle of a given radius with centre
 * (cx, cy)
 *
 * Returns 1 if (x-cx)^2 + (y-cy)^2 <= radius^2, 0 else
 *
 * @param x: x coordinate
 * @param y: y coordinate
 * @param cx: x coordinate of the centre of the circle
 * @param cy: y coordinate of the centre of the circle
 * @param radius: radius of the circle
 */
int in_circle(long x, long y, long cx, long cy, float radius) {
  return powf(x - cx, 2.0f) + powf(y - cy, 2.0f) <= radius * radius;
}

/**
 * \brief Check if the value is inside the bound.
 *
 * The value is inside the bound if i>=1 and i<=nx
 *
 * @param i: value to be checked
 * @param nx: maximum value
 */
int out_of_bounds(long i, long nx) { return i < 1 || i > nx; }

/**
 * \brief Check if the coordinate is inside the image bounds.
 *
 * @param i: x coordinate of the value to be checked
 * @param j: y coordinate of the value to be checked
 * @param nx: maximum value in x dimension
 * @param ny: maximum value in y dimension
 */
int in_image(long i, long j, long nx, long ny) {
  return !out_of_bounds(i, nx) && !out_of_bounds(j, ny);
}

/**
 * \brief Calculate the number of pixels inside a circle
 *
 * @param radius: radius of the circle
 */
int gauss_circle(float radius) {
  long l_max = ceil(radius) + 1;
  long n = 0;
  for (long i = -l_max; i <= l_max; ++i) {
    for (long j = -l_max; j <= l_max; ++j) {
      if (in_circle(i, j, 0, 0, radius)) {
        ++n;
      }
    }
  }
  return n;
}
