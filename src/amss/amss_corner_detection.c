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

/*
   Compute the differential operator L(u) as stated in Alvarez' paper
   */
void curvature(matrix u, /* original image: unchanged! */
    long nx,             /* image dimension in x direction */
    long ny,             /* image dimension in y direction */
    float h,             /* grid size; assumed quadratic! */
    matrix curv)         /* on exit, returns curvature values for each pixel */

{
    for (long i = 1; i <= nx; ++i) {
        for (long j = 1; j <= ny; ++j) {
        }
    }
    dummies(curv, nx, ny);
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
/*
       draws corners into the image u, at locations specified by v
       v[i][j]=255 for corners;
       v[i][j]=0   else;
       */

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

    return;

} /* draw_corners */

/*
   Find the solution to the linear system corresponding to the least squares approximation
   for each corner.
   */
void least_squares_approx(list_ptr chains, vector t, float q, long nx, long ny)
{
    matrix M; /* System matrix */
    vector b; /* Value vector */
    vector x; /* Solution of linear system */

    int n_iter = chains->chain_length;

    /* !!! CARE: Iteration in least squares starts with 1 !!! */
    alloc_matrix(&M, n_iter + 1, 3);
    alloc_vector(&b, n_iter + 1);
    alloc_vector(&x, 3);

    /* Initialise system matrix.
       System matrix stays the same across all corners */
    for (long k = 0; k < n_iter; ++k) {
        M[k + 1][1] = powf(t[k], 0.75f) - powf(t[0], 0.75f);
        M[k + 1][2] = 1;
    }

    for (node_ptr current = list_head(chains); current != NULL; current = current->next) {
        chain_ptr current_chain = current->chain;

        /* Initialise value vector */
        for (long k = 0; k < n_iter; ++k) {
            /* Distance to first point of corner chain */
            float diffx = current_chain[k].x - current_chain[0].x;
            float diffy = current_chain[k].y - current_chain[0].y;
            b[k + 1] = sqrt(diffx * diffx + diffy * diffy);
        }

        /* Approximate line using least squares approximation */
        solve_normal_equations(n_iter, 2, M, b, x);

        /* Need to correct slope from estimation as stated in the paper */
        current->slope = 0.071917f + 1.029484f * x[1];

        /* If slope is too small, remove chain because in this case, the corner evolution is
         * independent of time, which is not what we want. */
        if (fabs(current->slope) < 0.1) {
            list_delete(chains, current);
            continue;
        }

        current->angle = degrees(2.f * atan(1.f / powf(current->slope, 2.0f)));

        /* Calculate approximation error
         * FIXME: Is the error calculation correct? */
        float error = 0;
        for (long k = 1; k < n_iter + 1; ++k) {
            float acc = b[k] - (current->slope * M[k][1] + x[2]);
            error += acc * acc;
        }
        current->error = error / (n_iter + 1);

        /* Reconstruct initial corner position by extrapolation using
         * the corner bisector unit vector and slope */
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

    disalloc_matrix(M, n_iter + 1, 3);
    disalloc_vector(b, n_iter + 1);
    disalloc_vector(x, 3);
}

/*
   Detect corners in the given image using the AMSS.
*/
list_ptr amss_corner_detection

    (matrix u,      /* original image !! gets altered !! */
        long nx,    /* image dimension in x direction */
        long ny,    /* image dimension in y direction */
        float ht,   /* timestep size (<= 0.1) */
        long t_0,   /* initial scale for corner detection */
        long t_max, /* max scale */
        long step,  /* step size for scale space exploration */
        float q,    /* quantile of pixels to keep */
        matrix v)   /* output */

{
    vector t;         /* vector that contains the scales t0 ... tn of the scale space */
    matrix curv;      /* curvature of image at a certain timescale */
    matrix prev_amss; /* image at previous timescale */

    float h = 1.0;                       /* grid size; assuming we have a quadratic grid !*/
    float two_h_inv = 1.0f / (2.0f * h); /* Helper variable for faster computation */
    long n_iter = (t_max - t_0) / ht;

    printf("Number of iterations: %ld\n", n_iter);

    /* Linked list to store corner chains in */
    list_ptr chains = new_list(n_iter);

    alloc_vector(&t, n_iter);
    alloc_matrix(&curv, nx + 2, ny + 2);
    alloc_matrix(&prev_amss, nx + 2, ny + 2);

    t[0] = t_0;

    curvature(u, nx, ny, h, curv);
    write_pgm(curv, nx, ny, "init_curv.pgm", NULL);
    write_pgm(u, nx, ny, "init_amss.pgm", NULL);

    /* Use a higher threshold for the gradient in the first iteration */
    float kn = 0;
    while (kn < t_0) {
        amss(ht, nx, ny, h, h, u);
        kn += ht;
    }
    curvature(u, nx, ny, h, curv);

    write_pgm(curv, nx, ny, "t0_curv.pgm", NULL);
    write_pgm(u, nx, ny, "t0_amss.pgm", NULL);

    /* Find initial set of corners */
    for (long i = 1; i <= nx; i++) {
        for (long j = 1; j <= ny; j++) {
            /* if u[i][j] is bigger (or smaller) than all of its 8 neighbours
               then it is considered a local extremum */
            int max, min, all_zero;
            max = min = all_zero = 1;

            /* Check if curv[i][j] is a local extremum */
            for (long k = i - 1; k <= i + 1; ++k) {
                for (long l = j - 1; l <= j + 1; ++l) {
                    /* Can't be the maximum anymore */
                    if (curv[k][l] > curv[i][j]) {
                        max = 0;
                    }
                    /* Can't be the minimum anymore */
                    if (curv[k][l] < curv[i][j]) {
                        min = 0;
                    }
                    if (curv[k][l] != 0.0f) {
                        all_zero = 0;
                    }
                }
            }

            if ((min || max) && !all_zero) {
                node_ptr new = list_append_new(chains);
                new->chain[0] = new_pixel(i, j, curv[i][j]);
            }
        }
    }

    printf("Found %d extrema at the initial scale\n", list_size(chains));

    /* Track corners across time scales */
    for (long k = 1; k < n_iter; ++k) {
        printf("Sequences left: %d\n", chains->size);
        t[k] = t[0] + k * ht;

        imgcpy(u, prev_amss, nx + 2, ny + 2);

        /* First: calculate amss at time level t[k] in-place */
        amss(ht, nx, ny, h, h, u);

        char path[256];
        sprintf(path, "/home/daniel/Uni/Thesis/.tmp/amss%03ld.pgm", k);
        write_pgm(u, nx, ny, path, NULL);

        memset(path, '0', 256);
        sprintf(path, "/home/daniel/Uni/Thesis/.tmp/curv%03ld.pgm", k);
        write_pgm(curv, nx, ny, path, NULL);

        curvature(u, nx, ny, h, curv);

        /* Find new maximum for each chain at new timescale */
        for (node_ptr current = list_head(chains); current != NULL; current = current->next) {
            /* Position of corner at last scale */
            pixel_t last = current->chain[k - 1];
            long i = last.x;
            long j = last.y;

            assert(i >= 1);
            assert(j >= 1);

            /* Compute first order derivatives by means of central differences */
            float dx = (prev_amss[i - 1][j] - prev_amss[i + 1][j]) * two_h_inv;
            float dy = (prev_amss[i][j - 1] - prev_amss[i][j + 1]) * two_h_inv;

            /* Norm gradient to 1 */
            float norm_inv = 1 / (sqrt(dx * dx + dy * dy) + FLT_EPSILON);
            dx *= norm_inv;
            dy *= norm_inv;

            /* New locations for extrema search */
            long idx = (long)round(i + dx);
            long jdy = (long)round(j + dy);

            pixel_t cmax = new_pixel(-1, -1, 0);

            /* Search for new curvature extremum in 8-neighbourhood around pixel (idx, jdy)
               as a starting point for the next iteration */
            for (long l = idx - 1; l <= idx + 1; ++l) {
                for (long m = jdy - 1; m <= jdy + 1; ++m) {
                    /* We are not interested in these values since we can't track
                       any further in the next iteration if we're already at 0 */
                    if (l < 1 || l > nx || m < 1 || m > ny)
                        continue;
                    if (fabs(curv[l][m]) > fabs(cmax.curv)) {
                        cmax.x = l;
                        cmax.y = m;
                        cmax.curv = curv[l][m];
                    }
                }
            }

            /* Apply criteria to filter out spurious corner sequences
               and add non-spurious corners to the chain */
            int valid = cmax.x > -1 && cmax.y > -1;
            int same_sign = sgn(current->chain[k - 1].curv) == sgn(cmax.curv);
            int large_enough = sgn(cmax.curv) > CURV_EPSILON;

            if (valid && same_sign && large_enough) {
                current->chain[k] = cmax;
            } else {
                list_delete(chains, current);
            }
        }
    }

    /* Compute linear model */
    least_squares_approx(chains, t, q, nx, ny);
    print_list(chains);

    /* Only keep quantile of corner chains */
    printf("Chain size: %d\n", chains->size);
    //error_percentile(chains, 0.1);
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

    disalloc_matrix(curv, nx + 2, ny + 2);
    disalloc_matrix(prev_amss, nx + 2, ny + 2);
    return chains;
}

list_ptr test_detection(char* in, char* out)
{
    matrix u, f, overlaid;
    long nx, ny;
    read_pgm_and_allocate_memory(in, &nx, &ny, &u);

    alloc_matrix(&overlaid, nx + 2, ny + 2);
    alloc_matrix(&f, nx + 2, ny + 2);

    imgcpy(u, overlaid, nx + 2, ny + 2);

    list_ptr chains = amss_corner_detection(u, nx, ny, 0.1, 1, 10, 1, 0.05f, f);

    if (out != NULL) {
        draw_corners(overlaid, nx, ny, f);
        write_pgm(overlaid, nx, ny, out, NULL);
    }

    disalloc_matrix(f, nx + 2, ny + 2);
    disalloc_matrix(overlaid, nx + 2, ny + 2);
    disalloc_matrix(u, nx + 2, ny + 2);
    return chains;
}
