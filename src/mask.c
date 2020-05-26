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
            if (out_of_bounds(k, nx) || out_of_bounds(l, ny) ||
                !in_circle(k, l, i, j, radius)) {
              continue;
            }
            u[k][l] = 255.0;
          }
        }
      }
    }
  }
}
