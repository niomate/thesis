#include "amss.h"
#include "../utils.h"
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

void amss

    (float ht,     /* time step size, 0 < ht <= 0.1 */
        long nx,   /* image dimension in x direction */
        long ny,   /* image dimension in y direction */
        float hx,  /* pixel size in x direction */
        float hy,  /* pixel size in y direction */
        float** u) /* input: original image ;  output: smoothed */

/*
 affine morphological scale-space, explicit scheme
*/

{
    long i, j; /* loop variables */
    float** f; /* u at old time level */
    float h = hy;

    alloc_matrix(&f, nx + 2, ny + 2);

    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
            f[i][j] = u[i][j];

    dummies(f, nx, ny);

    /* loop */
    for (i = 1; i <= nx; i++) {
        for (j = 1; j <= ny; j++) {
            /*float dx = (f[i + 1][j] - f[i - 1][j]) / 2*h;*/
            /*float dy = (f[i][j + 1] - f[i][j - 1]) / 2*h;*/
            /*float dxx = (f[i + 1][j] - 2.0 * f[i][j] + f[i - 1][j]) / h*h;*/
            /*float dyy = (f[i][j + 1] - 2.0 * f[i][j] + f[i][j - 1]) / h*h;*/
            /*float dxy;*/
            /*if (dx * dy < 0.0)*/
            /*dxy = (f[i + 1][j + 1] - f[i][j + 1] - f[i + 1][j] + f[i][j] + f[i - 1][j - 1] -*/
            /*f[i][j - 1] - f[i - 1][j] + f[i][j]) /*/
            /*2*h*h;*/
            /*else*/
            /*dxy = (-f[i - 1][j + 1] + f[i][j + 1] + f[i + 1][j] - f[i][j] - f[i + 1][j - 1] +*/
            /*f[i][j - 1] + f[i - 1][j] - f[i][j]) /*/
            /*2*h*h;*/

            /*[> evolution <]*/
            /*float help = dx * dx * dyy + dy * dy * dxx - 2.0 * dx * dy * dxy;*/
            /*u[i][j] = f[i][j] + ht * sgn(help) * pow (fabs(help), 0.33333333f);*/

            float dx = 2.f * (f[i + 1][j] - f[i - 1][j]) + f[i + 1][j + 1] - f[i - 1][j + 1] + f[i + 1][j - 1] - f[i - 1][j - 1];

            float dy = 2.f * (f[i][j + 1] - f[i][j - 1]) + f[i + 1][j + 1] - f[i + 1][j - 1] + f[i - 1][j + 1] - f[i - 1][j - 1];

            float l0, l1, l2, l3, l4;
            if (fabs(dx) >= fabs(dy)) {
                l0 = 0.25f * (2.f * dx * dx + dy * dy - fabs(dx * dy));
            } else {
                l0 = 0.25f * (2.f * dy * dy + dx * dx - fabs(dx * dy));
            }

            l1 = 2.f * l0 - dy * dy;
            l2 = 2.f * l0 - dx * dx;
            l3 = -l0 + 0.5f * (dx * dy + dx * dx + dy * dy);
            l4 = -l0 + 0.5f * (-dx * dy + dx * dx + dy * dy);

            float help = -4.0f * l0 * f[i][j] + l1 * (f[i][j + 1] + f[i][j - 1]) + l2 * (f[i + 1][j] + f[i - 1][j]) + l3 * (f[i + 1][j - 1] + f[i - 1][j + 1]) + l4 * (f[i + 1][j + 1] + f[i - 1][j - 1]);

            u[i][j] = f[i][j] + ht * sgn(help) * powf(fabs(help), 0.33333333f);
        }
    }

    /* ---- free memory for f ---- */

    disalloc_matrix(f, nx + 2, ny + 2);

    return;

} /* amss */
