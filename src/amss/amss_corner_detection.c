#include "../utils.h"
#include "amss.h"
#include "chain.h"
#include "least_squares.h"
#include <assert.h>
#include <float.h>
#include <getopt.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef float** matrix;
typedef float* vector;

float CURV_EPSILON = 0.01;

int DEBUG = 1;

#define sqr(x) (x * x)

/*
   Compute the differential operator L(u) as stated in Alvarez' paper
   */
void curvature(matrix u, /* original image: unchanged! */
    long nx,             /* image dimension in x direction */
    long ny,             /* image dimension in y direction */
    float h,             /* grid size; assumed quadratic! */
    matrix Lu)           /* on exit, returns curvature values for each pixel */

{
    matrix f;
    alloc_matrix(&f, nx + 2, ny + 2);
    imgcpy(u, f, nx + 2, ny + 2);
    dummies(f, nx, ny);

    for (long i = 1; i <= nx; ++i) {
        for (long j = 1; j <= ny; ++j) {
            /*float dx = (u[i + 1][j] - u[i - 1][j]) / 2 * h;*/
            /*float dy = (u[i][j + 1] - u[i][j - 1]) / 2 * h;*/
            /*float dxx = (u[i + 1][j] - 2.0 * u[i][j] + u[i - 1][j]) / h *
             * h;*/
            /*float dyy = (u[i][j + 1] - 2.0 * u[i][j] + u[i][j - 1]) / h *
             * h;*/
            /*float dxy;*/
            /*dxy = (u[i + 1][j + 1] - u[i][j + 1] - u[i + 1][j] + u[i][j]*/
            /*+ u[i - 1][j - 1] - u[i][j - 1] - u[i - 1][j] + u[i][j])*/
            /*/ 2 * h * h;*/
            /*dxy += (-u[i - 1][j + 1] + u[i][j + 1] + u[i + 1][j] - u[i][j]*/
            /*- u[i + 1][j - 1] + u[i][j - 1] + u[i - 1][j] - u[i][j])*/
            /*/ 2 * h * h;*/
            /*dxy /= 2.0f;*/
            /*float help = dx * dx * dyy + dy * dy * dxx - 2.0 * dx * dy *
             * dxy;*/
            /*Lu[i][j] = help;*/
            float dx = 2.f * (f[i + 1][j] - f[i - 1][j]) + f[i + 1][j + 1] - f[i - 1][j + 1]
                + f[i + 1][j - 1] - f[i - 1][j - 1];

            float dy = 2.f * (f[i][j + 1] - f[i][j - 1]) + f[i + 1][j + 1] - f[i + 1][j - 1]
                + f[i - 1][j + 1] - f[i - 1][j - 1];

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

            float help = -4.0f * l0 * f[i][j] + l1 * (f[i][j + 1] + f[i][j - 1])
                + l2 * (f[i + 1][j] + f[i - 1][j]) + l3 * (f[i + 1][j - 1] + f[i - 1][j + 1])
                + l4 * (f[i + 1][j + 1] + f[i - 1][j - 1]);
            Lu[i][j] = help;
        }
    }
}

void write_normalized(matrix u, long nx, long ny, char* name)
{
    matrix f;
    alloc_matrix(&f, nx + 2, ny + 2);

    imgcpy(u, f, nx + 2, ny + 2);

    float max, min;
    max = -INFINITY;
    min = INFINITY;

    for (int i = 1; i <= nx; ++i) {
        for (int j = 1; j <= ny; ++j) {
            if (f[i][j] > max)
                max = f[i][j];
            if (f[i][j] < min)
                min = f[i][j];
        }
    }

    for (int i = 1; i <= nx; ++i) {
        for (int j = 1; j <= ny; ++j) {
            if (max > 255 && min < 0) {
                f[i][j] = 255.f * (f[i][j] - min) / (max - min);
            } else if (min >= 0) {
                f[i][j] = min + (f[i][j] - min) * (255.f - min) / (max - min);
            } else {
                f[i][j] = max * (f[i][j] - min) / (max - min);
            }
        }
    }

    write_pgm(f, nx, ny, name, NULL);
    disalloc_matrix(f, nx, ny);
}

void draw_corners(matrix u, /* image, altered ! */
    long nx,                /* image dimension in x direction */
    long ny,                /* image dimension in y direction */
    matrix v)               /* image with corner locations */
{
    long i, j; /* loop variables */
    long k, l; /* loop variables */

    for (i = 5; i <= nx - 4; i++)
        for (j = 5; j <= ny - 4; j++)
            if (v[i][j] == 255.0)
                /* draw corner */
                for (k = i - 4; k <= i + 4; k++)
                    for (l = j - 4; l <= j + 4; l++) {
                        /* black outer circle */
                        if ((k - i) * (k - i) + (l - j) * (l - j) <= 20)
                            u[k][l] = 0.0;
                        /* white interior */
                        if ((k - i) * (k - i) + (l - j) * (l - j) <= 6)
                            u[k][l] = 255.0;
                    }

} /* draw_corners */

/*
   Find the solution to the linear system corresponding to the least squares
   approximation for each corner.
*/
void least_squares_approx(list_ptr chains, vector t, float q, long nx, long ny)
{
    vector d;
    vector tau;

    int n_iter = chains->chain_length;

    alloc_vector(&d, n_iter);
    alloc_vector(&tau, n_iter);

    for (long k = 0; k < n_iter; ++k) {
        tau[k] = pow(t[k], 0.75f) - powf(t[0], 0.75f);
    }

    for (node_ptr current = list_head(chains); current != NULL; current = current->next) {
        chain_ptr current_chain = current->chain;

        /* Initialise value vector */
        for (long k = 0; k < n_iter; ++k) {
            /* Distance to first point of corner chain */
            float diffx = current_chain[k].x - current_chain[0].x;
            float diffy = current_chain[k].y - current_chain[0].y;
            d[k] = sqrt(sqr(diffx) + sqr(diffy));
        }

        float a, b, r_sqr;

        /* Approximate line using linear regression */
        if (!linear_regression(tau, d, &a, &b, &r_sqr, n_iter)) {
            list_delete(chains, current);
            continue;
        }

        current->slope = a;
        current->angle = degrees(2.f * atan(1.f / powf(current->slope, 2.0f)));
        current->error = r_sqr;

        /* Reconstruct initial corner position by extrapolation using
         * the corner bisector unit vector and slope */
        /* Review: Verify that this actually works! */
        float tdiff = pow(t[n_iter - 1], 0.75) - pow(t[0], 0.75);
        float ux = (current_chain[n_iter - 1].x - current_chain[0].x) / tdiff;
        float uy = (current_chain[n_iter - 1].y - current_chain[0].y) / tdiff;
        float tpow = pow(t[0], 0.75);
        float x0 = current_chain[0].x - ux * tpow;
        float y0 = current_chain[0].y - uy * tpow;
        current->corner_tip = new_pixel(x0, y0, NAN);
    }

    if (DEBUG) {
        printf("Found %d corner sequences.\n", list_size(chains));
    }

    disalloc_vector(tau, n_iter);
    disalloc_vector(d, n_iter);
}

/*
   Detect corners in the given image using the AMSS.
*/
list_ptr amss_corner_detection

    (matrix f,      /* original image !! gets altered !! */
        long nx,    /* image dimension in x direction */
        long ny,    /* image dimension in y direction */
        float ht,   /* timestep size (<= 0.1) */
        long t_0,   /* initial scale for corner detection */
        long t_max, /* max scale */
        float q,    /* quantile of pixels to keep */
        matrix v)   /* output */

{
    vector t; /* vector that contains the scales t0 ... tn of the scale space */
    matrix u;
    matrix Lu;     /* curvature of image at a certain timescale */
    matrix prev_u; /* image at previous timescale */

    float h = 1.0;                       /* grid size; assuming we have a quadratic grid !*/
    float two_h_inv = 1.0f / (2.0f * h); /* Helper variable for faster computation */
    long n_iter = (t_max - t_0) / ht;

    printf("Number of iterations: %ld\n", n_iter);

    /* Linked list to store corner chains in */
    list_ptr chains = new_list(n_iter);

    alloc_vector(&t, n_iter);
    alloc_matrix(&u, nx + 2, ny + 2);
    alloc_matrix(&Lu, nx + 2, ny + 2);
    alloc_matrix(&prev_u, nx + 2, ny + 2);

    dummies(f, nx, ny);

    imgcpy(f, u, nx + 2, ny + 2);

    t[0] = t_0;

    float kn = 0;
    while (kn < t_0) {
        amss(ht, nx, ny, h, h, u);
        kn += ht;
    }
    curvature(u, nx, ny, h, Lu);

    write_normalized(Lu, nx, ny, "initcurv.pgm");

    int max, min;
    /* Find initial set of corners */
    for (long i = 1; i <= nx; i++) {
        for (long j = 1; j <= ny; j++) {
            /* if u[i][j] is bigger (or smaller) than all of its 8 neighbours
               then it is considered a local extremum */
            max = 1;
            min = 1;
            for (int l = i - 1; l <= i + 1; ++l) {
                for (int m = j - 1; m <= j + 1; ++m) {
                    if (l == i && m == j)
                        continue;
                    max &= Lu[i][j] > Lu[l][m];
                    min &= Lu[i][j] < Lu[l][m];
                }
            }

            if (min || max) {
                node_ptr new = list_append_new(chains);
                new->chain[0] = new_pixel(i, j, Lu[i][j]);
            }
        }
    }

    /************************************** DEBUG ********************************************/
    matrix c;
    alloc_matrix(&c, nx + 2, ny + 2);
    for (int i = 1; i <= nx; ++i) {
        for (int j = 1; j <= ny; ++j) {
            c[i][j] = 0.0;
        }
    }
    printf("Found %d extrema at the initial scale\n", list_size(chains));
    for (node_ptr current = chains->head; current != NULL; current = current->next) {
        //printf("Extremum found at (%ld, %ld)\n", current->chain[0].x, current->chain[0].y);
        c[current->chain[0].x][current->chain[0].y] = 255.0;
    }
    write_pgm(c, nx, ny, "initcorners.pgm", NULL);
    disalloc_matrix(c, nx + 2, ny + 2);
    /************************************** DEBUG ********************************************/

    /* Track corners across time scales */
    for (long k = 1; k < n_iter; ++k) {
        //printf("\n############################[Iteration %ld]###########################\n\n", k);
        t[k] = t[0] + k * ht;

        imgcpy(u, prev_u, nx + 2, ny + 2);

        /* Since we already computed the differential operator, we just need
         * to perform the evolution on the previous timelevel */
        for (int i = 1; i <= nx; ++i) {
            for (int j = 1; j <= ny; ++j) {
                u[i][j] = u[i][j] + ht * sgn(Lu[i][j]) * pow(fabs(Lu[i][j]), 0.333333333);
            }
        }
        curvature(u, nx, ny, h, Lu);
        dummies(Lu, nx, ny);

        /************************************** DEBUG ********************************************/
        char path[256];
        sprintf(path, "/home/daniel/Uni/Thesis/.tmp/amss%03ld.pgm", k);
        write_pgm(u, nx, ny, path, NULL);

        memset(path, '0', 256);
        sprintf(path, "/home/daniel/Uni/Thesis/.tmp/curv%03ld.pgm", k);
        write_pgm(Lu, nx, ny, path, NULL);
        /************************************** DEBUG ********************************************/

        node_ptr current = list_head(chains);
        /* Find new maximum for each chain at new timescale */
        while (current != NULL) {
            /* Position of corner at last scale */
            long i = current->chain[k - 1].x;
            long j = current->chain[k - 1].y;

            assert(i >= 1);
            assert(j >= 1);

            /* Compute first order derivatives by means of central differences
             */
            float dx = (prev_u[i - 1][j] - prev_u[i + 1][j]) * two_h_inv;
            float dy = (prev_u[i][j - 1] - prev_u[i][j + 1]) * two_h_inv;

            /* Norm gradient to 1 */
            float norm_inv = 1 / (sqrt(dx * dx + dy * dy) + FLT_EPSILON);
            dx *= norm_inv;
            dy *= norm_inv;

            /* New locations for extrema search */
            /*long idx = (long)round(i + dx);*/
            /*long jdy = (long)round(j + dy);*/
            long idx = (long)(i + dx);
            long jdy = (long)(j + dy);

            pixel_t cmax = new_pixel(-1, -1, 0);

            /* Search for new curvature extremum in 8-neighbourhood around
               pixel (idx, jdy) as a starting point for the next iteration */
            for (long l = idx - 1; l <= idx + 1; ++l) {
                for (long m = jdy - 1; m <= jdy + 1; ++m) {
                    /* We are not interested in these values since we can't
                       track any further in the next iteration if we're already at 0 */
                    if (l < 1 || l > nx || m < 1 || m > ny)
                        continue;
                    if (fabs(Lu[l][m]) > fabs(cmax.curv)) {
                        cmax.x = l;
                        cmax.y = m;
                        cmax.curv = Lu[l][m];
                    }
                }
            }

            /* Apply criteria to filter out spurious corner sequences
               and add non-spurious corners to the chain */
            int valid = cmax.x > -1 && cmax.y > -1;
            int same_sign = sgn(current->chain[k - 1].curv) == sgn(cmax.curv);
            int large_enough = fabs(cmax.curv) > CURV_EPSILON;

            // print_pixel(cmax);
            // printf(", old curv: %f", current->chain[k - 1].curv);
            // printf(", %f\n", cmax.curv);
            // printf("Same sign %d, Large enough %d\n", same_sign, large_enough);
            if (valid && same_sign && large_enough) {
                current->chain[k] = cmax;
            } else {
                list_delete(chains, current);
            }
            current = current->next;
        }
    }

    /* Compute linear model */
    least_squares_approx(chains, t, q, nx, ny);

    print_list_minimal(chains);

    /* Only keep quantile of corner chains */
    printf("Chain size: %d\n", chains->size);
    error_percentile(chains, q);
    printf("Total corners: %d\n", chains->size);

    for (int i = 1; i < nx; ++i) {
        for (int j = 1; j < ny; ++j) {
            v[i][j] = 0.0f;
        }
    }

    /* Draw a marker at each corner position (only if valid) */
    printf("Corner locations: \n");
    for (node_ptr current = list_head(chains); current != NULL; current = current->next) {
        long x, y;
        x = current->corner_tip.x;
        y = current->corner_tip.y;
        printf("(%ld, %ld)\n", x, y);
        if (x < 0 || y < 0 || x > nx || y > ny) {
            continue;
        }
        v[x][y] = 255.0;
    }

    disalloc_matrix(Lu, nx + 2, ny + 2);
    disalloc_matrix(prev_u, nx + 2, ny + 2);
    return chains;
}

list_ptr test_detection(char* in, char* out)
{
    matrix f, v;
    long nx, ny;
    read_pgm_and_allocate_memory(in, &nx, &ny, &f);

    alloc_matrix(&v, nx + 2, ny + 2);

    list_ptr chains = amss_corner_detection(f, nx, ny, 0.1, 1, 20, 0.05, v);

    if (out != NULL) {
        draw_corners(f, nx, ny, v);
        write_pgm(f, nx, ny, out, NULL);
    }

    disalloc_matrix(f, nx + 2, ny + 2);
    disalloc_matrix(v, nx + 2, ny + 2);
    return chains;
}
