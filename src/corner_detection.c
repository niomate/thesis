#include "utils.h"
#include <stdarg.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define pi 3.1415927

#define MAX(x, y) ((x) > (y) ? (x) : (y))

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*               CORNER DETECTION WITH THE STRUCTURE TENSOR                 */
/*                                                                          */
/*                  (Copyright Joachim Weickert, 12/2005)                   */
/*                                                                          */
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                  Modified by Jana Gusenburger for                      */
/*                use in their bachelor's thesis titled                     */
/*      "Exploring circular corner regions in PDE-based inpainting"         */
/*                       (October '19 - June '20)                           */
/*--------------------------------------------------------------------------*/
/*                                                                          */
/* Modifications include:                                                   */
/*  - Command line interface                                                */
/*  - Adapted non maximum suppression                                       */
/*  - Dynamic thresholding for cornerness measure                           */
/*  - Mask computation                                                      */
/*                                                                          */
/*--------------------------------------------------------------------------*/

/**
 * Enum for convolution types for better readability */
typedef enum {
  GAUSSIAN = 0, /*!< Gaussian convolution */
  PILLBOX = 1   /*!< "Pillbox" kernel convolution-> simple averaging in disk
                   shaped neighbourhood */
} convolution;

/**
 * Enum for cornerness measure for better readability */
typedef enum { ROHR = 0, TOMASI = 1, HARRIS = 2 } corner_detector;
const char *CORNER_DETECTOR_NAMES[3] = {"Rohr", "Tomasi-Kanade",
                                        "Foerstner-Harris"};

/**
 * Enum for boundary conditions for better readability */
typedef enum {
  DIRICHLET = 0, /*!< Dirichlet boundary conditions,
                   i.e. 0 outside of the known grid*/
  NEUMANN = 1,   /*!< Neumann boundary conditions,
                   i.e. reflecting boundary conditions*/
  PERIOD = 2     /*!< Periodic boundary conditions*/
} boundary_condition;

/**
 * Additional options for corner detection extension */
typedef struct options {
  bool enable_cnms; /*<! Flag to enable circular non maximum suppression. If set
                       to false, just computes regular non maximum suppression*/
  bool enable_tppt; /*<! Flag to enable total pixel percentage thresholding. If
                       set to false, compute regular percentile thresholding on
                       cornerness value*/
  bool display_mask; /*<! Flag to enable computing the mask. If enabled, the
                        output of the programme will be an image that represents
                        the inpainting mask */
  convolution conv;  /*<! Convolution type */
} options_t;

options_t OPTIONS =
    (options_t){false, false, false, GAUSSIAN}; /*<! Default options */

/**
 * Struct containing information about a pixel. Contains position and value.
 * Used in total_pixel_percentage_thresholding to get a more accurate result.
 */
typedef struct {
  long i;
  long j;
  float val;
} pixel_t;

/**
 * \brief Float compare function.
 */
int float_cmp(const void *a, const void *b) {
  const float v1 = *((const float *)a);
  const float v2 = *((const float *)b);
  if (v1 < v2)
    return -1;
  return v1 > v2;
}

/**
 * \brief Pixel compare function.
 */
int pixel_cmp(const void *a, const void *b) {
  const pixel_t v1 = *((const pixel_t *)a);
  const pixel_t v2 = *((const pixel_t *)b);
  if (v1.val < v2.val)
    return -1;
  return v1.val > v2.val;
}


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

/**
 * \brief Uniform average in disk shaped neighbourhood of size radius.
 *
 * Uses Dirichlet boundary conditions.
 * Convolve the image with a normalised disk shaped kernel for a sharper
 * representation of the influence area of the non-vanishing part of a Gaussian
 * kernel.
 *
 * @param radius; size of the kernel
 * @param nx: size of the image in x direction
 * @param ny: size of the image in y direction
 * @param hx: grid size of the image in x direction
 * @param hy: grid size of the image in y direction
 */
void pillbox_conv(float radius, long nx, long ny, float hx, float hy,
                  float **f) {
  long i, j, k, l; /* Loop variables */
  float **help;    /* Copy of the image */
  alloc_matrix(&help, nx + 2, ny + 2);

  long length = ceil(radius) + 1; /* Length of the kernel in either direction */
  long lmax = 2 * length + 1;     /* Total length of the kernel in positive and
                                     negative diurection combined */

  float **kernel; /* Convolution kernel */
  alloc_matrix(&kernel, lmax, lmax);

  long n = 0; /* Number of pixels that are 1, used for normalization */

  /* Compute pillbox kernel for the given radius */
  for (i = 0; i < lmax; ++i) {
    for (j = 0; j < lmax; ++j) {
      /* Check if (i, j) is in a circle with origin
       * (length, length) */
      if (in_circle(i, j, length, length, radius)) {
        kernel[i][j] = 1;
        n++;
      } else {
        kernel[i][j] = 0;
      }
    }
  }

  float norm = 1.0 / (float)n;

  /* Normalise kernel */
  for (i = 0; i < lmax; ++i) {
    for (j = 0; j < lmax; ++j) {
      kernel[i][j] = kernel[i][j] * norm;
    }
  }

  /* Convolution step */
  for (i = 1; i <= nx; ++i) {
    for (j = 1; j <= ny; j++) {
      help[i][j] = 0;
      for (k = -length; k <= length; ++k) {
        for (l = -length; l <= length; ++l) {
          if (in_image(i + k, j + l, nx, ny))
            help[i][j] += f[i + k][j + l] * kernel[length + k][length + l];
        }
      }
    }
  }

  /* Put computed data back */
  imgcpy(help, f, nx + 2, ny + 2);

  disalloc_matrix(help, nx + 2, ny + 2);
  disalloc_matrix(kernel, lmax, lmax);
}

void gauss_conv

    (float sigma,           /* standard deviation of Gaussian */
     long nx,               /* image dimension in x direction */
     long ny,               /* image dimension in y direction */
     float hx,              /* pixel size in x direction */
     float hy,              /* pixel size in y direction */
     float precision,       /* cutoff at precision * sigma */
     boundary_condition bc, /* type of boundary condition */
                            /* 0=Dirichlet, 1=reflecing, 2=periodic */
     float **f)             /* input: original image ;  output: smoothed */

/*
 Gaussian convolution.
*/

{
  long i, j, p; /* loop variables */
  long length;  /* convolution vector: 0..length */
  float sum;    /* for summing up */
  float *conv;  /* convolution vector */
  float *help;  /* row or column with dummy boundaries */

  /* ------------------------ diffusion in x direction -------------------- */

  /* calculate length of convolution vector */
  length = (long)(precision * sigma / hx) + 1;
  if ((bc != 0) && (length > nx)) {
    printf("gauss_conv: sigma too large \n");
    exit(0);
  }

  /* allocate storage for convolution vector */
  alloc_vector(&conv, length + 1);

  /* calculate entries of convolution vector */
  for (i = 0; i <= length; i++)
    conv[i] = 1 / (sigma * sqrt(2.0 * 3.1415926)) *
              exp(-(i * i * hx * hx) / (2.0 * sigma * sigma));

  /* normalization */
  sum = conv[0];
  for (i = 1; i <= length; i++) {
    sum = sum + 2.0 * conv[i];
  }

  for (i = 0; i <= length; i++)
    conv[i] = conv[i] / sum;

  /* allocate storage for a row */
  alloc_vector(&help, nx + length + length);

  for (j = 1; j <= ny; j++) {
    /* copy in row vector */
    for (i = 1; i <= nx; i++)
      help[i + length - 1] = f[i][j];

    /* assign boundary conditions */
    if (bc == DIRICHLET) /* Dirichlet boundary conditions */
      for (p = 1; p <= length; p++) {
        help[length - p] = 0.0;
        help[nx + length - 1 + p] = 0.0;
      }
    else if (bc == NEUMANN) /* reflecting b.c. */
      for (p = 1; p <= length; p++) {
        help[length - p] = help[length + p - 1];
        help[nx + length - 1 + p] = help[nx + length - p];
      }
    else if (bc == PERIOD) /* periodic b.c. */
      for (p = 1; p <= length; p++) {
        help[length - p] = help[nx + length - p];
        help[nx + length - 1 + p] = help[length + p - 1];
      }

    /* convolution step */
    for (i = length; i <= nx + length - 1; i++) {
      /* calculate convolution */
      sum = conv[0] * help[i];
      for (p = 1; p <= length; p++)
        sum = sum + conv[p] * (help[i + p] + help[i - p]);
      /* write back */
      f[i - length + 1][j] = sum;
    }
  } /* for j */

  /* disallocate storage for a row */
  disalloc_vector(help, nx + length + length);

  /* disallocate convolution vector */
  disalloc_vector(conv, length + 1);

  /* ------------------------ diffusion in y direction -------------------- */

  /* calculate length of convolution vector */
  length = (long)(precision * sigma / hy) + 1;
  if ((bc != 0) && (length > ny)) {
    printf("gauss_conv: sigma too large \n");
    exit(0);
  }

  /* allocate storage for convolution vector */
  alloc_vector(&conv, length + 1);

  /* calculate entries of convolution vector */
  for (j = 0; j <= length; j++)
    conv[j] = 1 / (sigma * sqrt(2.0 * 3.1415927)) *
              exp(-(j * j * hy * hy) / (2.0 * sigma * sigma));

  /* normalization */
  sum = conv[0];
  for (j = 1; j <= length; j++)
    sum = sum + 2.0 * conv[j];
  for (j = 0; j <= length; j++)
    conv[j] = conv[j] / sum;

  /* allocate storage for a row */
  alloc_vector(&help, ny + length + length);

  for (i = 1; i <= nx; i++) {
    /* copy in column vector */
    for (j = 1; j <= ny; j++)
      help[j + length - 1] = f[i][j];

    /* assign boundary conditions */
    if (bc == DIRICHLET) /* Dirichlet boundary conditions */
      for (p = 1; p <= length; p++) {
        help[length - p] = 0.0;
        help[ny + length - 1 + p] = 0.0;
      }
    else if (bc == NEUMANN) /* reflecting b.c. */
      for (p = 1; p <= length; p++) {
        help[length - p] = help[length + p - 1];
        help[ny + length - 1 + p] = help[ny + length - p];
      }
    else if (bc == PERIOD) /* periodic b.c. */
      for (p = 1; p <= length; p++) {
        help[length - p] = help[ny + length - p];
        help[ny + length - 1 + p] = help[length + p - 1];
      }

    /* convolution step */
    for (j = length; j <= ny + length - 1; j++) {
      /* calculate convolution */
      sum = conv[0] * help[j];
      for (p = 1; p <= length; p++)
        sum = sum + conv[p] * (help[j + p] + help[j - p]);
      /* write back */
      f[i][j - length + 1] = sum;
    }
  } /* for i */

  /* disallocate storage for a row */
  disalloc_vector(help, ny + length + length);

  /* disallocate convolution vector */
  disalloc_vector(conv, length + 1);

  return;

} /* gauss_conv */

void struct_tensor

    (float **v,   /* image !! gets smoothed on exit !! */
     long nx,     /* image dimension in x direction */
     long ny,     /* image dimension in y direction */
     float hx,    /* pixel size in x direction */
     float hy,    /* pixel size in y direction */
     float sigma, /* noise scale */
     float rho,   /* integration scale */
     float **dxx, /* element of structure tensor, output */
     float **dxy, /* element of structure tensor, output */
     float **dyy) /* element of structure tensor, output */

/*
 Calculates the structure tensor.
*/

{
  long i, j;            /* loop variables */
  float dv_dx, dv_dy;   /* derivatives of v */
  float h1, h2, h3, h4; /* time saver */

  /* ---- initialization ---- */

  for (i = 1; i <= nx; i++) {
    for (j = 1; j <= ny; j++) {
      dxx[i][j] = 0.0;
      dxy[i][j] = 0.0;
      dyy[i][j] = 0.0;
    }
  }

  float alpha = 0.2;
  h1 = (1.0 - 2.0 * alpha) / (2.0 * hx);
  h2 = alpha / (2.0 * hx);
  h3 = (1.0 - 2.0 * alpha) / (2.0 * hy);
  h4 = alpha / (2.0 * hy);

  /* ---- smoothing at noise scale ---- */

  if (sigma > 0.0) {
    gauss_conv(sigma, nx, ny, hx, hy, 3.0, 1, v); /* refl. b.c. */
  }

  /* ---- building tensor product ---- */

  dummies(v, nx, ny);

  for (i = 1; i <= nx; i++) {
    for (j = 1; j <= ny; j++) {
      /* modified Sobel derivatives */
      dv_dx = h1 * (v[i + 1][j] - v[i - 1][j]) +
              h2 * (v[i + 1][j - 1] - v[i - 1][j - 1] + v[i + 1][j + 1] -
                    v[i - 1][j + 1]);
      dv_dy = h3 * (v[i][j + 1] - v[i][j - 1]) +
              h4 * (v[i - 1][j + 1] - v[i - 1][j - 1] + v[i + 1][j + 1] -
                    v[i + 1][j - 1]);
      /* tensor product */
      dxx[i][j] = dv_dx * dv_dx;
      dxy[i][j] = dv_dx * dv_dy;
      dyy[i][j] = dv_dy * dv_dy;
    }
  }

  /* ---- smoothing at integration scale, Dirichlet b.c. ---- */

  if (rho > 0.0) {
    if (OPTIONS.conv == PILLBOX) {
      pillbox_conv(rho, nx, ny, hx, hy, dxx);
      pillbox_conv(rho, nx, ny, hx, hy, dxy);
      pillbox_conv(rho, nx, ny, hx, hy, dyy);
    } else if (OPTIONS.conv == GAUSSIAN) {
      gauss_conv(rho, nx, ny, hx, hy, 3.0, 0, dxx);
      gauss_conv(rho, nx, ny, hx, hy, 3.0, 0, dxy);
      gauss_conv(rho, nx, ny, hx, hy, 3.0, 0, dyy);
    }
  }

  return;

} /* struct_tensor */

void PA_trans

    (float a11,   /* coeffs of (2*2)-matrix */
     float a12,   /* coeffs of (2*2)-matrix */
     float a22,   /* coeffs of (2*2)-matrix */
     float *c,    /* 1. comp. of 1. eigenvector, output */
     float *s,    /* 2. comp. of 1. eigenvector, output */
     float *lam1, /* larger  eigenvalue, output */
     float *lam2) /* smaller eigenvalue, output */

/*
 Principal axis transformation, checked for correctness.
*/

{
  float help, norm; /* time savers */

  /* ---- compute eigenvalues and eigenvectors ---- */

  help = sqrt(powf(a22 - a11, 2.0) + 4.0 * a12 * a12);

  if (help == 0.0)
  /* isotropic situation, eigenvectors arbitrary */
  {
    *lam1 = *lam2 = a11;
    *c = 1.0;
    *s = 0.0;
  }

  else if (a11 > a22) {
    *lam1 = 0.5 * (a11 + a22 + help);
    *lam2 = 0.5 * (a11 + a22 - help);
    *c = a11 - a22 + help;
    *s = 2.0 * a12;
  }

  else {
    *lam1 = 0.5 * (a11 + a22 + help);
    *lam2 = 0.5 * (a11 + a22 - help);
    *c = 2.0 * a12;
    *s = a22 - a11 + help;
  }

  /* ---- normalize eigenvectors ---- */

  norm = sqrt(*c * *c + *s * *s);
  if (norm > 0.0) {
    *c = *c / norm;
    *s = *s / norm;
  } else {
    *c = 1.0;
    *s = 0.0;
  }

  return;

} /* PA_trans */

/**
 * \brief Computes the inpainting mask for the given image.
 *
 * Given the corner map v and a radius, draw the circular corner regions onto
 * the image u.
 *
 * @param u: the original image, altered afterwards
 * @param nx: size of the image in x direction
 * @param ny: size of the image in x direction
 * @param radius: radius of the corner regions
 * @param v: corner map
 *
 */
void mask(float **u, long nx, long ny, int radius, float **v) {
  long i, j, k, l; /* loop variables */

  int r = radius + 2; /* To create a circle, we have to start the iteration
                         outside of it to not get a weird looking rectangle */

  for (i = 1; i <= nx; ++i) {
    for (j = 1; j <= ny; ++j) {
      u[i][j] = 0.0f;
    }
  }

  for (i = 1; i <= nx; ++i) {
    for (j = 1; j <= ny; j++) {
      if (v[i][j] != 255.0) {
        /* Not a corner */
        continue;
      } else {
        /* Create a circle of radius r around the current pixel */
        for (k = i - r; k <= i + r; k++) {
          for (l = j - r; l <= j + r; l++) {
            if (out_of_bounds(k, nx) || out_of_bounds(l, ny) ||
                !in_circle(k, l, i, j, radius)) {
              continue;
            }
            u[k][l] = 255.0f;
          }
        }
      }
    }
  }
}

/**
 * \brief Search for the local maxima of the cornerness measure.
 *
 * Given a map of cornerness values, find the local maxima inside an
 * 8-neighbourhood and discard every value that does not comply. This gives a
 * less noisy and more accurate approximation of the actual positions of the
 * corners.
 *
 * @param w cornerness map, unchanged
 * @param nx: size of the image in x direction
 * @param ny: size of the image in x direction
 * @param v: corner map, afterwards contains the values of corners that are
 * being kept, 0 elsewhere
 */
void non_maximum_suppression(float **w, long nx, long ny, float **v) {
  long i, j, k, l; /* Loop variables */
  for (i = 1; i <= nx; i++) {
    for (j = 1; j <= ny; j++) {
      if (w[i][j] <= 0.001) {
        continue;
      }
      int max = 1;
      /* Check if current pixel is maximum in 3x3 neighbourhood */
      for (k = i - 1; k <= i + 1; ++k) {
        for (l = j - 1; l <= j + 1; ++l) {
          if (w[k][l] > w[i][j]) {
            max = 0;
          }
        }
      }
      if (max) {
        v[i][j] = w[i][j];
      } else {
        v[i][j] = 0.0f;
      }
    }
  }
}

/**
 * \brief Discard corners that are already covered by 'better' corners.
 *
 * Iterate over all cornerness values and set the values of corners that are not
 * the maximum inside a disk shaped neighbourdhood of the given radius to 0.
 * The idea is that combined with the total pixel percentage thresholding, we
 * can distribute the centres of the mask pixels more evenly across the image by
 * thinning out the amount of corners in denser areas by discarding those
 * corners that are already inside the disk of another corner.
 * This implementation is not focused on speed or efficiency, it only serves as
 * a proof of concept to show that an approach like this makes sense.
 *
 * This implementation also contains some bugs in edge case scenarios.
 * e.g.
 * if a corner is found to be covered by another corner, it is discarded. If
 * however this next corner is also found to be covered by an even better
 * corner, it is also discarded and we completely lose the information of the
 * first corner.
 *
 * Currently, I don't know how to fix this edge case. If this proves to be a
 * successful concept, then this would be one of the first points of
 * improvement.
 *
 * @param w cornerness map, unchanged
 * @param nx: size of the image in x direction
 * @param ny: size of the image in x direction
 * @param radius; mask size, i.e. radius of the disk surrounding each corner
 * @param v: corner map, afterwards contains the values of corners that are
 * being kept, 0 elsewhere
 */
void non_maximum_suppression_circle(float **w, long nx, long ny, long radius,
                                    float **v) {
  long i, j, k, l; /* Loop variables */
  int max;         /* Variables for non maximum suppression */

  for (i = 1; i <= nx; i++) {
    for (j = 1; j <= ny; j++) {
      max = 1;
      /* To speed up the process, weed out outliers by first comparing them to
       * an epsilon */
      if (w[i][j] <= 0.001) {
        continue;
      }
      /* Check if current pixel is maximum in circular neighbourhood of given
       * radius */
      for (k = i - radius; k <= i + radius; k++) {
        for (l = j - radius; l <= j + radius; l++) {
          if (!in_image(k, l, nx, ny) || !in_circle(k, l, i, j, radius)) {
            continue;
          }
          if (w[k][l] > w[i][j]) {
            max = 0;
          }
        }
      }
      if (max) {
        v[i][j] = w[i][j];
      } else {
        v[i][j] = 0;
      }
    }
  }
}

/**
 * \brief Using a cornerness map, find the corners that should be kept
 * in order to obtain a mask containing a certain percentage of pixels.
 *
 * Compute how many corners with a mask disk of size radius can be kept such
 * that the resulting mask contains at most the given percentage of pixels.
 * This method is more accurate in computing masks with a similar amount
 * of pixels than using the percentile approach.
 * Furthermore, it prevents the case where the all corners have the same value
 * and thus choosing a percentile too low removes all corners from the map.
 * In this approach we sort the corners by their cornerness value and then, as
 * long as we have space left in the mask, put them into the mask one by one.
 * This gives a better approximation to the desired mask size.
 *
 * @param u: cornerness map, thresholded afterwards, 255 at corner locations, 0
 *           elsewhere
 * @param nx: size of the image in x direction
 * @param ny: size of the image in x direction
 * @param radius; mask size, i.e. radius of the disk surrounding each corner
 * @param perc: percentage of corners that should be kept
 */
void total_pixel_percentage_thresholding(float **u, long nx, long ny,
                                         int radius, float perc) {
  long i, j;
  long n_total = nx * ny;
  pixel_t vals[n_total];

  /* Compute the area explicitly to get the most accurate result. */
  float max_area_p_corner = gauss_circle(radius);

  long n_corners = floor((perc * n_total) / max_area_p_corner);
  printf("Area in pixels per corner: %f\n", max_area_p_corner);
  printf("Maximum number of corners allowed: %ld\n", n_corners);
  printf("Maximum percentage possible: %f\n",
         100.0f * (max_area_p_corner * n_corners) / n_total);

  /* Flatten cornerness map and prepare for sorting */
  long idx = 0;
  for (i = 1; i <= nx; ++i) {
    for (j = 1; j <= ny; ++j) {
      if (u[i][j] > 0.01) {
        vals[idx++] = (pixel_t){i, j, u[i][j]};
      }
    }
  }

  qsort(vals, idx, sizeof(pixel_t), pixel_cmp);

  /* Initialize */
  for (i = 1; i <= nx; i++) {
    for (j = 1; j <= ny; j++) {
      u[i][j] = 0;
    }
  }

  /* Keep best corners */
  long nums = 0;
  for (i = idx - 1; i >= MAX(idx - n_corners, 0); --i) {
    long x = vals[i].i;
    long y = vals[i].j;
    u[x][y] = 255.0f;
    nums++;
  }

  printf("Inserted %ld corners\n", nums);
}

/**
 * \brief Computes the n-th percentile given by the parameter perc on a
 * cornerness map
 *
 * Compute the n-th percentile on the given cornerness map and threshold it
 * against the threshold computed this way. Afterwards, u is 0 at locations
 * that were below the threshold and 255 for locations above the threshold.
 *
 * The parameter perc specifies how many percent of corners should be **kept**.
 *
 * @param u: cornerness map, thresholded afterwards, 255 at corner locations, 0
 *           elsewhere
 * @param nx: size of the image in x direction
 * @param ny: size of the image in x direction
 * @param perc: percentage of corners that should be kept
 */
void percentile_thresholding(float **u, long nx, long ny, float perc) {
  long i, j;
  long n_total = nx * ny;
  pixel_t vals[n_total];

  /* Flatten cornerness map and prepare for sorting */
  long idx = 0;
  for (i = 1; i <= nx; ++i) {
    for (j = 1; j <= ny; ++j) {
      if (u[i][j] > 0.01) {
        vals[idx++] = (pixel_t){i, j, u[i][j]};
      }
    }
  }

  long index = clamp((long)ceil(idx * perc), 0, idx);
  qsort(vals, idx, sizeof(pixel_t), pixel_cmp);

  /* Initialize */
  for (i = 1; i <= nx; i++) {
    for (j = 1; j <= ny; j++) {
      u[i][j] = 0;
    }
  }

  /* Keep best corners */
  for (i = idx - 1; i >= MAX(idx - index, 0); --i) {
    long x = vals[i].i;
    long y = vals[i].j;
    u[x][y] = 255.0f;
  }
}

/**
 * \brief Corner detection using the structure tensor
 *
 * Find image features using a structure tensor based corner detection method
 * like Harris, Rohr or Tomasi-Kanade.
 * Incorporates some extensions like modified non maximum suppression and
 * percentile thresholding.
 *
 * General order of actions is:
 *  - Presmooth the image
 *  - Compute structure tensor information
 *  - Apply the cornerness measure
 *  - Apply non maximum suppression
 *  - Apply percentile thresholding or TPPT
 *
 * @param f original image, unaltered
 * @param v corner locations, 255 at locations of image feature, 0 elsewhere
 * @param nx size of image in x direction
 * @param ny size of image in y direction
 * @param hx grid size of image in x direction
 * @param hy grid size of image in y direction
 * @param type type of cornerness measure, a.k.a corner detector
 * @param perc percentile parameter for percentile thresholding
 * @param radius radius of the disk in the mask image and size of disk for CNMS
 * @param sigma noise scale, i.e. scale at which the image is smoothed before
 *        running the corner detection
 * @param rho integration scale, i.e. size of the Gaussian kernel that is used
 *        to average directional information from the structure tensor
 */
void corner_detection(float **f, float **v, long nx, long ny, float hx,
                      float hy, corner_detector type, float perc, float radius,
                      float sigma, float rho) {
  long i, j;                 /* loop variables */
  float **u;                 /* loop */
  float **dxx, **dxy, **dyy; /* tensor components */
  float **w;                 /* field for corner detection */
  float c, s;                /* cosine, sine */
  float lam1, lam2;          /* eigenvalues */
  float det, trace;          /* determinant, trace */

  /* allocate storage */
  alloc_matrix(&u, nx + 2, ny + 2);
  alloc_matrix(&dxx, nx + 2, ny + 2);
  alloc_matrix(&dxy, nx + 2, ny + 2);
  alloc_matrix(&dyy, nx + 2, ny + 2);
  alloc_matrix(&w, nx + 2, ny + 2);

  /* copy f into u */
  imgcpy(f, u, nx + 2, ny + 2);

  /* calculate structure tensor of u */
  struct_tensor(u, nx, ny, hx, hy, sigma, rho, dxx, dxy, dyy);

  if (type == ROHR)
    /* Rohr: det(J) */
    for (i = 1; i <= nx; i++) {
      for (j = 1; j <= ny; j++) {
        w[i][j] = dxx[i][j] * dyy[i][j] - dxy[i][j] * dxy[i][j];
      }
    }

  if (type == TOMASI)
    /* Tomasi-Kanade: smaller eigenvalue */
    for (i = 1; i <= nx; i++) {
      for (j = 1; j <= ny; j++) {
        PA_trans(dxx[i][j], dxy[i][j], dyy[i][j], &c, &s, &lam1, &lam2);
        w[i][j] = lam2;
      }
    }

  if (type == HARRIS)
    /* Foerstner-Harris: det(J)/trace(J) */
    for (i = 1; i <= nx; i++) {
      for (j = 1; j <= ny; j++) {
        /* Adding a small epsilon to prevent division by zero */
        trace = dxx[i][j] + dyy[i][j] + 0.000001;
        det = dxx[i][j] * dyy[i][j] - dxy[i][j] * dxy[i][j];
        w[i][j] = 2 * det / trace;
      }
    }

  dummies(w, nx, ny);

  if (!OPTIONS.enable_cnms) {
    non_maximum_suppression(w, nx, ny, v);
  } else {
    non_maximum_suppression_circle(w, nx, ny, radius, v);
  }

  if (!OPTIONS.enable_tppt) {
    percentile_thresholding(v, nx, ny, perc);
  } else {
    total_pixel_percentage_thresholding(v, nx, ny, radius, perc);
  }

  /* free storage */
  disalloc_matrix(u, nx + 2, ny + 2);
  disalloc_matrix(dxx, nx + 2, ny + 2);
  disalloc_matrix(dxy, nx + 2, ny + 2);
  disalloc_matrix(dyy, nx + 2, ny + 2);
  disalloc_matrix(w, nx + 2, ny + 2);
} /* corner_detection */

/*--------------------------------------------------------------------------*/

/**
 * \brief Highlight corner locations in the original image by drawing a cross at
 * the corner location.
 *
 * @param u original image, highlighted corner locations afterwards
 * @param nx size of image in x direction
 * @param ny size of image in y direction
 * @param v corner locations
 */
void draw_corners_x(float **u, long nx, long ny, float **v) {
  long i, j; /* loop variables */
  long k, l; /* loop variables */

  dummies(u, nx, ny);

  for (i = 5; i <= nx - 4; i++) {
    for (j = 5; j <= ny - 4; j++) {
      /* draw corner */
      if (v[i][j] == 255.0) {
        for (k = -5; k <= 5; k++) {
          if (k == 0) {
            u[i][j] = 255 - u[i][j];
          } else {
            u[i + k][j] = 255 - u[i + k][j];
            u[i][j + k] = 255 - u[i][j + k];
          }
        }
      }
    }
  }
} /* draw_corners */

/**
 * \brief Highlight corner locations in the original image by drawing a black
 * and white circle at the corner location.
 *
 * @param u original image, highlighted corner locations afterwards
 * @param nx size of image in x direction
 * @param ny size of image in y direction
 * @param v corner locations
 */
void draw_corners_o(float **u, long nx, long ny, float **v) {
  long i, j; /* loop variables */
  long k, l; /* loop variables */

  dummies(u, nx, ny);

  for (i = 5; i <= nx - 4; i++) {
    for (j = 5; j <= ny - 4; j++) {
      if (v[i][j] == 255.0) {
        /* draw corner */
        for (k = i - 4; k <= i + 4; k++) {
          for (l = j - 4; l <= j + 4; l++) {
            /* black outer circle */
            if ((k - i) * (k - i) + (l - j) * (l - j) <= 20)
              u[k][l] = 0.0;
            /* white interior */
            if ((k - i) * (k - i) + (l - j) * (l - j) <= 6)
              u[k][l] = 255.0;
          }
        }
      }
    }
  }
} /* draw_corners */

int main(int argc, char **argv) {
  char *in;    /* Filename input image */
  float **u;   /* Input image */
  float **v;   /* Image with corner location */
  long nx, ny; /* Image size in x, y direction */
  long count;  /* Number of corners */

  /* Variables given from command line interface */
  char *out = "corners.pgm";     /* Filename output image */
  corner_detector type = HARRIS; /* Type of corner detector */
  float q = 0.1;                 /* Percentile parameter to compute threshold */
  float sigma = 1;               /* Noise scale */
  float rho = 2.5;               /* Integration scale */
  float mask_factor = 0; /* Factor to scale integration scale with to obtain
                            radius of corner regions */
  char display_type = 'o';

  printf("\n");
  printf("CORNER DETECTION WITH THE STRUCTURE TENSOR\n\n");

  int c;
  while ((c = getopt(argc, argv, "q:c:s:r:o:m:d:CDMP")) != -1) {
    switch (c) {
    case 'q':
      q = atof(optarg);
      if (q < 0 || q > 1) {
        fprintf(stderr, "Percentile has to be between 0 and 1!\n");
        abort();
      }
      break;
    case 'c':
      type = (corner_detector)atoi(optarg);
      break;
    case 's':
      sigma = atof(optarg);
      break;
    case 'r':
      rho = atof(optarg);
      break;
    case 'o':
      out = optarg;
      break;
    case 'm':
      mask_factor = atof(optarg);
      break;
    case 'C':
      OPTIONS.enable_cnms = true;
      break;
    case 'P':
      OPTIONS.enable_tppt = true;
      break;
    case 'D':
      OPTIONS.conv = PILLBOX;
      break;
    case 'M':
      OPTIONS.display_mask = true;
      break;
    case 'd':
      printf("Arg: %c\n", *optarg);
      if (*optarg != 'x' && *optarg != 'o') {
        fprintf(stderr, "Invalid corner display type! Defaulting to circle\n");
      } else {
        display_type = *optarg;
      }
      printf("Using %c to mark corner locations.\n", display_type);
      break;
    default:
      abort();
    }
  }

  if (optind == argc) {
    fprintf(stderr, "No input file specified. Aborting...\n");
    abort();
  }

  if (argc - optind > 1) {
    fprintf(stderr, "Too many positional arguments. Aborting...\n");
    abort();
  }

  in = argv[optind];

  read_pgm_and_allocate_memory(in, &nx, &ny, &u);

  /* Allocate storage for corner locations */
  alloc_matrix(&v, nx + 2, ny + 2);

  /* Initialize corner map */
  for (long i = 1; i <= nx; ++i) {
    for (long j = 1; j <= ny; ++j) {
      v[i][j] = 0.0f;
    }
  }

  /* Compute mask radius in dependence of rho */
  float mask_radius = mask_factor; // mask_factor * rho;

  printf("Using parameters:\n");
  printf("  Corner detector:    %-s\n", CORNER_DETECTOR_NAMES[type]);
  printf("  Initial image:      %-s\n", in);
  printf("  Output image:       %-s\n", out);
  printf("  Sigma:              %-8.2f\n", sigma);
  printf("  Rho:                %-8.2f\n", rho);
  printf("  Percentile:         %-8.2f\n", q);
  printf("  Mask factor:        %-8.2f\n", mask_factor);
  printf("  Mask radius:        %-8.2f\n", mask_radius);
  printf("  CNMS enabled:       %-5s\n",
         OPTIONS.enable_cnms ? "true" : "false");
  printf("  TPPT enabled:       %-5s\n",
         OPTIONS.enable_tppt ? "true" : "false");
  printf("  Compute mask:       %-5s\n",
         OPTIONS.display_mask ? "true" : "false");
  printf("\n");

  corner_detection(u, v, nx, ny, 1.0, 1.0, type, q, mask_radius, sigma, rho);

  long n_corners = 0;
  for (long i = 1; i <= nx; ++i) {
    for (long j = 1; j <= ny; ++j) {
      if (v[i][j] == 255.0f) {
        n_corners++;
      }
    }
  }

  printf("Number of detected corners:       %-8ld\n", n_corners);

  if (OPTIONS.display_mask) {
    /* Create mask with disk shaped corner regions */
    mask(u, nx, ny, mask_radius, v);

    long n_mask_pixels = 0;
    for (long i = 1; i <= nx; ++i) {
      for (long j = 1; j <= ny; ++j) {
        if (u[i][j] == 255.0f) {
          n_mask_pixels++;
        }
      }
    }
    printf("Number of pixels in mask:       %-8ld\n", n_mask_pixels);
    printf("Percentage:                     %-8.2f\n",
           100.0f * ((float)n_mask_pixels) / (nx * ny));
  } else {
    /* Insert corner marks into the image */
    if (display_type == 'x')
      draw_corners_x(u, nx, ny, v);
    else
      draw_corners_o(u, nx, ny, v);
  }

  write_pgm(u, nx, ny, out, NULL);

  /* ---- disallocate storage ---- */

  disalloc_matrix(u, nx + 2, ny + 2);
  disalloc_matrix(v, nx + 2, ny + 2);
}
