#include "utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define pi 3.1415927

void mask

    (float **u,  /* image, altered ! */
     long nx,    /* image dimension in x direction */
     long ny,    /* image dimension in y direction */
     int radius, /* Radius of mask */
     float **v)  /* image with corner locations */

{
  long i, j; /* loop variables */
  long k, l; /* loop variables */
  dummies(u, nx, ny);

  int r = 2 * radius;
  for (i = 1; i <= nx; ++i) {
    for (j = 1; j <= ny; j++) {
      u[i][j] = 0;
    }
  }

  for (i = 1; i <= nx; ++i) {
    for (j = 1; j <= ny; j++) {
      if (v[i][j] != 255.0) {
        continue;
      }
      for (k = i - r; k <= i + r; k++) {
        for (l = j - r; l <= j + r; l++) {
          if (k < 1 || k > nx || l < 1 || l > ny)
            continue;
          float dist = powf(k - i, 2) + powf(l - j, 2);
          if (dist > radius * radius) {
            continue;
          }
          u[k][l] = 255.0;
        }
      }
    }
  }
}

/* Generate a random float between 0 and max */
float randomf(float max) { return ((float)rand() / (float)(RAND_MAX)) * max; }

void sample_mask

    (float **u, long nx, long ny, int radius, float perc, float **v)

/* Same as mask, but sample disk to keep perc% of pixels inside each disk
 * instead of a solid disk
 * */

{

  long i, j;
  long n;
  int total = radius * pi * pi;

  dummies(u, nx, ny);

  for (i = 1; i <= nx; ++i) {
    for (j = 1; j <= ny; ++j) {
      u[i][j] = 0;
    }
  }

  for (i = 1; i <= nx; ++i) {
    for (j = 1; j <= ny; ++j) {
      if (v[i][j] != 255.0) {
        continue;
      }
      /* Sample disk uniformly */
      for (n = 0; n < floor(perc * total); ++n) {
        /* Choose a random radius and angle */
        float r = randomf(radius);
        float phi = randomf(2 * pi);
        int x = i + round(r * cos(phi));
        int y = j + round(r * sin(phi));
        if (x < 1 || x > nx || y < 1 || y > ny) {
          continue;
        }
        u[x][y] = 255.0f;
      }
    }
  }
}
