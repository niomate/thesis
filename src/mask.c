#include "utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define pi 3.1415927

void create_initial_image_inpainting(float **u, long nx, long ny, float **mask) {
  for (long i = 0; i <= nx; ++i) {
    for (long j = 0; i <= ny; ++i) {
      if (mask[i][j] <= 0.0001) {
        u[i][j] = 0.0;
      }
    }
  }
}

void mask

    (float **u,  /* image, altered ! */
     long nx,    /* image dimension in x direction */
     long ny,    /* image dimension in y direction */
     int radius, /* Radius of mask */
     float **v)  /* image with corner locations */

{
  long i, j; /* loop variables */
  long k, l; /* loop variables */
  long n_mask_pixels = 0;
  dummies(u, nx, ny);

  int r = radius + 2; /* To create a circle, we have to start the iteration
                         outside of it to not get a weird looking rectangle */

  for (i = 1; i <= nx; ++i) {
    for (j = 1; j <= ny; j++) {
      u[i][j] = 0;
    }
  }

  for (i = 1; i <= nx; ++i) {
    for (j = 1; j <= ny; j++) {
      if (v[i][j] != 255.0) {
        /* Not a corner */
        continue;
      } else if (radius == 0) {
        /* Add only the current pixel to the mask since we do not use corner
         * regions */
        u[i][j] = 255.0;
      } else {
        /* Create a circle around of radius r around the current pixel */
        for (k = i - r; k <= i + r; k++) {
          for (l = j - r; l <= j + r; l++) {
            if (k < 1 || k > nx || l < 1 || l > ny)
              continue;
            float dist = powf(k - i, 2) + powf(l - j, 2);
            if (dist >= radius * radius) {
              continue;
            }
            u[k][l] = 255.0;
            n_mask_pixels++;
          }
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
