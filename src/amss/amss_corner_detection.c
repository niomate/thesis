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

typedef float **matrix;
typedef float *vector;

float CURV_EPSILON = 0.01;

int DEBUG = 1;

/*
    Compute the differential operator L(u) as stated in Alvarez' paper
*/
void curvature

(matrix u,    /* original image: unchanged! */
 long nx,     /* image dimension in x direction */
 long ny,     /* image dimension in y direction */
 float h,     /* grid size; assumed quadratic! */
 matrix curv) /* on exit, returns curvature values for each pixel */

{
    long i, j;

    matrix f; /* copy of original image */

    float two_h = 2.0 * h;
    float h_sqr = h * h;
    float two_h_sqr = 2.0 * h * h;

    alloc_matrix (&f, nx + 2, ny + 2);
    for (i = 1; i <= nx; ++i) {
        for (j = 1; j <= ny; j++) {
            f[i][j] = u[i][j];
        }
    }
    dummies (f, nx, ny);

    for (i = 1; i <= nx; ++i) {
        for (j = 1; j <= ny; ++j) {
            // Compute differential operator as proposed in Modelli et Ciomaga 2011
            float dx = 2.f * (u[i + 1][j] - u[i - 1][j]) + u[i + 1][j + 1] - u[i - 1][j + 1] +
                       u[i + 1][j - 1] - u[i - 1][j - 1];

            float dy = 2.f * (u[i][j + 1] - u[i][j - 1]) + u[i + 1][j + 1] - u[i + 1][j - 1] +
                       u[i - 1][j + 1] - u[i - 1][j - 1];

            float l0, l1, l2, l3, l4;
            if (fabs (dx) >= fabs (dy)) {
                l0 = 0.25f * (2.f * dx * dx + dy * dy - fabs (dx * dy));
            } else {
                l0 = 0.25f * (2.f * dy * dy + dx * dx - fabs (dx * dy));
            }

            l1 = 2.f * l0 - dy * dy;
            l2 = 2.f * l0 - dx * dx;
            l3 = -l0 + 0.5f * (dx * dy + dx * dx + dy * dy);
            l4 = -l0 + 0.5f * (-dx * dy + dx * dx + dy * dy);

            float help = -4.0f * l0 * u[i][j] + l1 * (u[i][j + 1] + u[i][j - 1]) +
                         l2 * (u[i + 1][j] + u[i - 1][j]) + l3 * (u[i + 1][j - 1] + u[i - 1][j + 1]) +
                         l4 * (u[i + 1][j + 1] + u[i - 1][j - 1]);

            curv[i][j] = sgn (help) * powf (fabs (help), 0.33333333f);
        }
    }
    dummies (curv, nx, ny);

    disalloc_matrix (f, nx + 2, ny + 2);
}


void iter_amss (matrix u, float kmax, long nx, long ny, long h, float ht) {
    float kn = 0;
    while (kn < kmax) {
        amss (ht, nx, ny, h, h, u);
        kn += ht;
    }
}

void write_normalized (matrix u, long nx, long ny, char *name) {
    matrix f;
    alloc_matrix (&f, nx + 2, ny + 2);

    imgcpy (u, f, nx + 2, ny + 2);

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
            f[i][j] = (f[i][j] - min) / (max - min) * 255.f;
        }
    }

    write_pgm (f, nx, ny, name, NULL);
    disalloc_matrix (f, nx, ny);
}


/*
    Find the solution to the linear system corresponding to the least squares approximation
    for each corner.
*/
void least_squares_approx (list_ptr chains, vector t, float q, long nx, long ny, matrix corners) {
    matrix M; /* System matrix */
    vector b; /* Value vector */
    vector x; /* Solution of linear system */

    int n_iter = chains->chain_length;

    /* !!! CARE: Iteration in least squares starts with 1 !!! */
    alloc_matrix (&M, n_iter + 1, 3);
    alloc_vector (&b, n_iter + 1);
    alloc_vector (&x, 3);

    /* Initialise system matrix.
    System matrix stays the same across all corners */
    for (long k = 0; k < n_iter; ++k) {
        M[k + 1][1] = powf (t[k], 0.75f) - powf (t[0], 0.75f);
        M[k + 1][2] = 1;
    }

    for (node_ptr current = list_head (chains); current != NULL; current = current->next) {
        chain_ptr current_chain = current->chain;

        /* Initialise value vector */
        for (long k = 0; k < n_iter; ++k) {
            /* Distance to first point of corner chain */
            float diffx = current_chain[k].x - current_chain[0].x;
            float diffy = current_chain[k].y - current_chain[0].y;
            b[k + 1] = sqrt (diffx * diffx + diffy * diffy);
        }

        /* Approximate line using least squares approximation */
        solve_normal_equations (n_iter, 2, M, b, x);

        /* Need to correct slope from estimation as stated in the paper */
        current->slope = 0.071917f + 1.029484f * x[1];

        /* If slope is too small, remove chain because in this case, the corner evolution is
         * independent of time, which is not what we want. */
        if (fabs (current->slope) < 0.1) {
            list_delete (chains, current);
            continue;
        }

        // printf ("Current slope (%ld, %ld): %f\n", current_chain[0].x, current_chain[0].y, current->slope);
        current->angle = degrees (2.f * atan (1.f / powf (current->slope, 2.0f)));

        // if (DEBUG)
        //     printf ("\nSlope: %f, Angle: %f\n", current->slope, current->angle);

        /* Calculate approximation error
         * FIXME: Is the error calculation correct? */
        float error = 0;
        for (long k = 1; k < n_iter + 1; ++k) {
            float acc = b[k] - (current->slope * M[k][1] + x[2]);
            error += acc * acc;
        }
        current->error = error / (n_iter + 1);

        // if (DEBUG)
        // printf ("Current error: %f\n", current->error);

        float diffx = (current_chain[0].x - current_chain[n_iter - 1].x) / b[n_iter];
        float diffy = (current_chain[0].y - current_chain[n_iter - 1].y) / b[n_iter];

        // printf ("Corner bisector (%ld, %ld): (%f, %f)\n", current_chain[0].x, current_chain[0].y, diffx, diffy);

        /* Reconstruct initial corner position by extrapolation using
        the corner bisector unit vector and slope */
        float helper = current->slope * powf (1.33333333f * t[0], 0.75f);

        // printf ("%f\n", helper);
        // printf ("%f, %f\n", diffx, diffy);

        current->corner_tip = new_pixel (round (current_chain[0].x - helper * diffx),
                                         round (current_chain[0].y - helper * diffy), NAN);

        // printf ("Computed tip: (%ld, %ld)\n\n", current->corner_tip.x, current->corner_tip.y);
    }

    if (DEBUG) {
        printf ("Found %d corner sequences.\n", list_size (chains));
    }

    disalloc_matrix (M, n_iter + 1, 3);
    disalloc_vector (b, n_iter + 1);
    disalloc_vector (x, 3);
}

/*
    Detect corners in the given image using the AMSS.
*/
list_ptr amss_corner_detection

(matrix u,   /* original image !! gets altered !! */
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

    printf ("Number of iterations: %ld\n", n_iter);

    /* Linked list to store corner chains in */
    list_ptr chains = new_list (n_iter);

    alloc_vector (&t, n_iter);
    alloc_matrix (&curv, nx + 2, ny + 2);
    alloc_matrix (&prev_amss, nx + 2, ny + 2);

    t[0] = t_0;

    curvature (u, nx, ny, h, curv);
    write_pgm (curv, nx, ny, "init_curv.pgm", NULL);
    write_pgm (u, nx, ny, "init_amss.pgm", NULL);


    /* Use a higher threshold for the gradient in the first iteration */
    iter_amss (u, t_0, nx, ny, h, ht);
    curvature (u, nx, ny, h, curv);

    write_pgm (curv, nx, ny, "t0_curv.pgm", NULL);
    write_pgm (u, nx, ny, "t0_amss.pgm", NULL);


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
                node_ptr new = list_append_new (chains);
                new->chain[0] = new_pixel (i, j, curv[i][j]);
            }
        }
    }

    printf ("Found %d extrema at the initial scale\n", list_size (chains));

    /* Track corners across time scales */
    for (long k = 1; k < n_iter; ++k) {
        t[k] = t[0] + k * ht;

        /* FIXME: This might not be necessary, could compute it in the previous step and store it in the pixel struct */
        imgcpy (u, prev_amss, nx + 2, ny + 2);

        /* First: calculate amss at time level t[k] in-place */
        amss (ht, nx, ny, h, h, u);

        char path[256];
        sprintf (path, "/home/danielg/uni/thesis/.tmp/amss%03ld.pgm", k);
        write_pgm (u, nx, ny, path, NULL);

        memset (path, '0', 256);
        sprintf (path, "/home/danielg/uni/thesis/.tmp/curv%03ld.pgm", k);
        write_pgm (curv, nx, ny, path, NULL);

        curvature (u, nx, ny, h, curv);

        /* Find new maximum for each chain at new timescale */
        for (node_ptr current = list_head (chains); current != NULL; current = current->next) {
            /* Position of corner at last scale */
            pixel_t last = current->chain[k - 1];
            long i = last.x;
            long j = last.y;

            assert (i >= 1);
            assert (j >= 1);

            /* Compute first order derivatives by means of central differences */
            float dx = (prev_amss[i - 1][j] - prev_amss[i + 1][j]) * two_h_inv;
            float dy = (prev_amss[i][j - 1] - prev_amss[i][j + 1]) * two_h_inv;

            /* Norm gradient to 1 */
            float norm_inv = 1 / (sqrt (dx * dx + dy * dy) + FLT_EPSILON);
            dx *= norm_inv;
            dy *= norm_inv;

            /* New locations for extrema search */
            long idx = (long)round (i + dx);
            long jdy = (long)round (j + dy);

            pixel_t cmax = new_pixel (-1, -1, 0);

            /* Search for new curvature extremum in 8-neighbourhood around pixel (idx, jdy)
            as a starting point for the next iteration */
            for (long l = idx - 1; l <= idx + 1; ++l) {
                for (long m = jdy - 1; m <= jdy + 1; ++m) {
                    /* We are not interested in these values since we can't track
                    any further in the next iteration if we're already at 0 */
                    if (l < 1 || l > nx || m < 1 || m > ny)
                        continue;
                    if (fabs (curv[l][m]) > fabs (cmax.curv)) {
                        cmax.x = l;
                        cmax.y = m;
                        cmax.curv = curv[l][m];
                    }
                }
            }

            /* Apply criteria to filter out spurious corner sequences
            and add non-spurious corners to the chain */
            int valid = cmax.x > -1 && cmax.y > -1;
            int same_sign = sgn (current->chain[k - 1].curv) == sgn (cmax.curv);
            int large_enough = sgn (cmax.curv) > CURV_EPSILON;

            if (valid && same_sign && large_enough) {
                current->chain[k] = cmax;
            } else {
                list_delete (chains, current);
            }
        }
    }
   //print_list(chains);

    // write_pgm (u, nx, ny, "tmax_amss.pgm", NULL);
    // write_pgm (curv, nx, ny, "tmax_curv.pgm", NULL);
    /* Compute linear model */
    least_squares_approx (chains, t, q, nx, ny, v);

    /* Only keep quantile of corner chains */
    error_percentile (chains, 0.1);

    //print_list (chains);

    for (long i = 1; i <= nx; ++i) {
        for (long j = 1; j <= ny; ++j) {
            v[i][j] = 0;
        }
    }

    for (node_ptr current = list_head (chains); current != NULL; current = current->next) {
        long x = current->corner_tip.x;
        long y = current->corner_tip.y;

        printf("Found corner at (%ld, %ld)\n", x, y);
        v[x][y] = 255.0;
    }


    if (DEBUG) {
        printf ("Corners left after quantile application: %d .\n", list_size (chains));
        /* TODO: visualize the corner angles usign Lemma 1 from the paper */
        write_pgm (v, nx, ny, "/home/danielg/uni/thesis/.tmp/corners.pgm", NULL);

        /* Print chains to image for better visualisation*/
        matrix chain_vis;
        alloc_matrix (&chain_vis, nx + 2, ny + 2);

        for (long i = 0; i < nx; ++i) {
            for (long j = 0; j < ny; ++j) {
                chain_vis[i][j] = 0;
            }
        }

        node_ptr current = list_head (chains);
        while (current != NULL) {
            chain_ptr c = current->chain;
            for (long k = 0; k < n_iter; ++k) {
                chain_vis[c[k].x][c[k].y] = 255;
            }
            current = current->next;
        }

        write_pgm (chain_vis, nx, ny, "/home/danielg/uni/thesis/.tmp/chains.pgm", NULL);
    }

    disalloc_vector (t, n_iter);
    disalloc_matrix (curv, nx + 2, ny + 2);
    disalloc_matrix (prev_amss, nx + 2, ny + 2);

    return chains;
}


list_ptr test_detection (char *in, char *out) {
    matrix u, f;
    long nx, ny;
    read_pgm_and_allocate_memory (in, &nx, &ny, &u);
    alloc_matrix (&f, nx + 2, ny + 2);
    list_ptr chains = amss_corner_detection (u, nx, ny, 0.1, 1, 20, 1, 0.05f, f);

    if (out != NULL)
        write_pgm (f, nx, ny, out, NULL);

    disalloc_matrix (f, nx + 2, ny + 2);
    disalloc_matrix (u, nx + 2, ny + 2);
    return chains;
}
