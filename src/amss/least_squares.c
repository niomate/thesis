#include "least_squares.h"
#include "../utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void pivot

(long n,    /* number of equations and unknowns, input */
 long j,    /* column for which we seach the pivot element */
 float **b, /* system matrix of size n * n, changed */
 float *d)  /* right hand side of size n, changed */

/*
 column pivot strategy:
 - searches the pivot element p in column j:
   |b[p][j]| >= |b[i][j]| for i=j,...,n
 - if j<>p, exchange the equations j and p
*/

{
    long i, k; /* loop variables */
    long p;    /* index of pivot element */
    float aux; /* auxiliary variable */


    /* ---- search pivot element ---- */

    p = j;
    for (i = j + 1; i <= n; i++)
        if (fabs (b[p][j]) >= fabs (b[i][j]))
            p = i;


    /* ---- exchange equations j and p  ---- */

    if (j != p) {
        /* exchange rows b[j][*] and b[p][*] */
        for (k = j; k <= n; k++) {
            aux = b[p][k];
            b[p][k] = b[j][k];
            b[j][k] = aux;
        }

        /* exchange right hand side entries d[j] and d[p] */
        aux = d[p];
        d[p] = d[j];
        d[j] = aux;
    }

    return;

} /* pivot */

/*--------------------------------------------------------------------------*/

void gauss

(long n,    /* number of equations, input */
 float **b, /* system matrix of size n * n, input */
 float *d,  /* right hand side of size n, input */
 float *x)  /* solution vector, output */

/*
 Gauss algorithm for solving a linear system of equations with pivoting
*/

{
    long i, j, k; /* loop variables */
    float sum;    /* for summing up */
    float factor; /* time saver */


    /* ---- bring linear system of size n * n in echelon form ---- */

    /* make sure that |b(1,1)| >= |b(k,1)| for k=2,...,n */
    pivot (n, 1, b, d);

    for (i = 1; i <= n - 1; i++) {
        /* make sure that |b(i,i)| >= |b(k,i)| for k=i+1,...,n */
        pivot (n, i, b, d);

        for (k = i + 1; k <= n; k++) {
            /* subtract b[k][i] / b[i][i] times equation i from equation k */
            factor = b[k][i] / b[i][i];
            for (j = i; j <= n; j++)
                b[k][j] = b[k][j] - factor * b[i][j];
            d[k] = d[k] - factor * d[i];
        }
    }

    /* ---- backward substitution ---- */

    for (i = n; i >= 1; i--) {
        sum = d[i];
        for (j = i + 1; j <= n; j++)
            sum = sum - b[i][j] * x[j];
        if (b[i][i] != 0.0)
            x[i] = sum / b[i][i];
        else
            printf ("singular linear system of equations\n");
    }

    return;

} /* gauss */

/*--------------------------------------------------------------------------*/

void solve_normal_equations

(long imax,  /* number of equations, input */
 long jmax,  /* number of unknowns, input */
 float **a,  /* system matrix of size imax * jmax, input */
 float *rhs, /* right hand side of size imax, input */
 float *w,   /* weight of each equation, input */
 float *x)   /* solution of normal equations, output */

/*
 solves the normal equations that arise from a linear system of imax
 equations with jmax unknowns, where equation i has the weight w[i]
*/

{
    long i, j, k; /* loop variables */
    float **b;    /* system matrix of normal system */
    float *d;     /* right hand side of normal system */


    /* ---- allocate storage ---- */

    alloc_matrix (&b, jmax + 1, jmax + 1);
    alloc_vector (&d, jmax + 1);


    /* ---- construct normal equations B x = d ---- */

    /* d = A^T W rhs */
    for (i = 1; i <= jmax; i++) {
        d[i] = 0.0;
        for (k = 1; k <= imax; k++)
            d[i] = d[i] + a[k][i] * w[k] * rhs[k];
    }

    /* B = A^T W A */
    for (i = 1; i <= jmax; i++)
        for (j = 1; j <= jmax; j++) {
            b[i][j] = 0.0;
            for (k = 1; k <= imax; k++)
                b[i][j] = b[i][j] + a[k][i] * w[k] * a[k][j];
        }


    /* ---- solve normal equations with Gauss algorithm ---- */

    gauss (jmax, b, d, x);


    /* ---- disallocate storage ---- */

    disalloc_matrix (b, jmax + 1, jmax + 1);
    disalloc_vector (d, jmax + 1);
    return;

} /* solve_normal_equations */

/*--------------------------------------------------------------------------*/
