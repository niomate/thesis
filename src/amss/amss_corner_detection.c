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

extern int DEBUG;

/*
    Compute the differential operator L(u) as stated in Alvarez' paper
*/
void curvature

(float **u,    /* original image: unchanged! */
 long nx,      /* image dimension in x direction */
 long ny,      /* image dimension in y direction */
 float h,      /* grid size; assumed quadratic! */
 float **curv) /* on exit, returns curvature values for each pixel */

{
    long i, j;
    float dx, dy, dxx, dxy, dyy;

    float **f; /* copy of original image */

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
            /* central spatial derivatives */
            dx = (f[i + 1][j] - f[i - 1][j]) / two_h;
            dy = (f[i][j + 1] - f[i][j - 1]) / two_h;
            dxx = (f[i + 1][j] - 2.0 * f[i][j] + f[i - 1][j]) / h_sqr;
            dyy = (f[i][j + 1] - 2.0 * f[i][j] + f[i][j - 1]) / h_sqr;

            dxy = (f[i + 1][j + 1] - f[i][j + 1] - f[i + 1][j] + f[i][j] + f[i - 1][j - 1] -
                   f[i][j - 1] - f[i - 1][j] + f[i][j]) /
                  two_h_sqr;
            dxy += (-f[i - 1][j + 1] + f[i][j + 1] + f[i + 1][j] - f[i][j] - f[i + 1][j - 1] +
                    f[i][j - 1] + f[i - 1][j] - f[i][j]) /
                   two_h_sqr;

            dxy /= 2.0;

            /* compute differential operator as proposed in Alvarez and Morales, 1997 */
            curv[i][j] = dx * dx * dyy + dy * dy * dxx - 2.0 * dx * dy * dxy;
        }
    }
}

void normalize_range

(float **u, /* image */
 long nx,   /* image dimension in x direction */
 long ny)   /* image dimension in y direction */

{
    long i, j; /* loop variables */
    float max, min;
    max = -INFINITY;
    min = INFINITY;

    for (i = 1; i <= nx; ++i) {
        for (j = 1; j <= ny; ++j) {
            if (u[i][j] > max) {
                max = u[i][j];
            }
            if (u[i][j] < min) {
                min = u[i][j];
            }
        }
    }

    for (i = 1; i <= nx; ++i) {
        for (j = 1; j <= ny; ++j) {
            u[i][j] = (u[i][j] - min) * 255.0f / (max - min);
        }
    }
}


void iter_amss (float **u, float target, long nx, long ny, long h, float ht) {
    float kn = 0;
    while (kn < target) {
        amss (ht, nx, ny, h, h, u);
        kn += ht;
    }
}


/*
    Find the solution to the linear system corresponding to the least squares approximation
    for each corner.
*/
void least_squares_approx (list_ptr chains, float *t, float q, long nx, long ny, float **corners) {
    float **M; /* System matrix */
    float *b;  /* Value vector */
    float *x;  /* Solution of linear system */

    int n_iter = list_chain_length (chains);

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
        if (DEBUG) {
            printf ("Current chain: ");
            print_chain (current_chain, n_iter);
        }

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
        current->angle = 2.f * atan (1.f / (current->slope * current->slope));

        if (DEBUG)
            printf ("\nSlope: %f, angle: %f\n", current->slope, current->angle);

        /* Calculate approximation error */
        float error = 0;
        for (long k = 1; k < n_iter + 1; ++k) {
            float acc = b[k] - (current->slope * M[k][1] + x[2]);
            error += acc * acc;
        }
        current->error = error / n_iter;

        if (DEBUG)
            printf ("Current error: %f\n\n", current->error);

        float diffx = current_chain[n_iter - 1].x - current_chain[0].x;
        float diffy = current_chain[n_iter - 1].y - current_chain[0].y;

        /* Reconstruct initial corner position by extrapolation using
        the corner bisector unit vector and slope */
        float norm_inv = 1.f / (b[n_iter] + FLT_EPSILON); /* b[n_iter] == sqrt(diffx²-diffy²) */
        float helper = current->slope * powf (1.33333333f * t[0], 0.75f) * norm_inv;

        current->corner_tip = pixel_new (round (current_chain[0].x - helper * diffx),
                                         round (current_chain[0].y - helper * diffy), NAN);
    }

    if (DEBUG) {
        printf ("Found %d corner sequences.\n", list_size (chains));
        printf ("Applying quantile %f to corner chains...\n", q);
    }

    /* Only keep quantile of corner chains */
    // threshold_error_quantile (chains, n_iter, n_corners_before, q);

    for (long i = 1; i <= nx; ++i) {
        for (long j = 1; j <= ny; ++j) {
            corners[i][j] = 0;
        }
    }

    for (node_ptr current = list_head (chains); current != NULL; current = node_next (current)) {
        long x = current->corner_tip.x;
        long y = current->corner_tip.y;
        corners[x][y] = 255.0;
    }

    if (DEBUG) {
        printf ("Corners left after quantile application: %d .\n", list_size (chains));
    }

    disalloc_matrix (M, n_iter + 1, 3);
    disalloc_vector (b, n_iter + 1);
    disalloc_vector (x, 3);
}

/*
    Detect corners in the given image using the AMSS.
*/
list_ptr amss_corner_detection

(float **u,  /* original image !! gets altered !! */
 long nx,    /* image dimension in x direction */
 long ny,    /* image dimension in y direction */
 float ht,   /* timestep size (<= 0.1) */
 long t_0,   /* initial scale for corner detection */
 long t_max, /* max scale */
 long step,  /* step size for scale space exploration */
 float q,    /* quantile of pixels to keep */
 float **v)  /* output */


{
    float *t;     /* vector that contains the scales t0 ... tn of the scale space */
    float **curv; /* curvature of image at a certain timescale */

    float h = 1.0;                       /* grid size; assuming we have a quadratic grid !*/
    float two_h_inv = 1.0f / (2.0f * h); /* Helper variable for faster computation */
    long n_iter = t_max / step + 1;      /* number of total iterations */

    /* Linked list to store corner chains in */
    list_ptr chains = list_new (n_iter);

    if (DEBUG)
        printf ("Number of iterations: %ld\n", n_iter);

    alloc_vector (&t, n_iter);
    alloc_matrix (&curv, nx + 2, ny + 2);

    if (DEBUG)
        printf ("Initial scale: %ld\n", t_0);

    t[0] = t_0;

    iter_amss (u, t_0, nx, ny, h, ht);

    if (DEBUG)
        write_pgm (u, nx, ny, "/home/danielg/uni/thesis/.tmp/curv_init.pgm", NULL);

    curvature (u, nx, ny, h, curv);
    dummies (curv, nx, ny);

    /* Find initial set of corners */
    for (long i = 1; i <= nx; i++) {
        for (long j = 1; j <= ny; j++) {
            /* if u[i][j] is bigger (or smaller) than all of its 8 neighbours
            then it is considered a local extremum (currently only comparing it to
            the 4 neighbours up, down, left and right */
            int extremum = 1;
            float zero_acc = 1;
            /* FIXME: This is super ugly */
            for (long k = -1; k <= 1; ++k) {
                for (long l = -1; l <= 1; ++l) {
                    if (!k && !l) // k == l == 0
                        continue;
                    if (fabs (curv[i][j]) < fabs (curv[i + k][j + l])) {
                        extremum = 0;
                    }
                    zero_acc *= curv[i + k][j + l];
                }
            }
            if (extremum && zero_acc != 0) {
                node_ptr new = list_insert (chains, NULL);
                node_add_to_chain (new, 0, i, j, curv[i][j]);
            }
        }
    }

    if (DEBUG)
        printf ("Found %d extrema at the initial scale\n", list_size (chains));

    /* Track corners across time scales */
    for (long k = 1; k < n_iter; ++k) {
        if (DEBUG)
            printf ("Iteration %ld\n", k);
        t[k] = t[k - 1] + step;

        /* First: calculate amss at time level t[k] in-place */
        iter_amss (u, step, nx, ny, h, ht);
        curvature (u, nx, ny, h, curv);
        dummies (curv, nx, ny);

        /* Find new maximum for each chain at new timescale */
        for (node_ptr current = list_head (chains); current != NULL; current = node_next (current)) {
            /* Position of corner at last scale */
            pixel_t last = current->chain[k - 1];
            long i = last.x;
            long j = last.y;

            assert (i >= 1);
            assert (j >= 1);

            /* Compute first order derivatives by means of central differences */
            float dx = (u[i - 1][j] - u[i + 1][j]) * two_h_inv;
            float dy = (u[i][j - 1] - u[i][j + 1]) * two_h_inv;

            /* Norm gradient to 1 */
            float norm_inv = 1 / (sqrt (dx * dx + dy * dy) + FLT_EPSILON);
            dx *= norm_inv;
            dy *= norm_inv;

            /* New locations for extrema search */
            long idx = (long)round (i + dx);
            long jdy = (long)round (j + dy);

            long imax = -1;
            long jmax = -1;
            long cmax = -1;

            /* Might have to have a look and see if this is actually correct */
            /* Search for new curvature extremum in 8-neighbourhood around pixel (idx, jdy) */
            /* FIXME: This is super ugly as well */
            for (long l = clamp (idx - 1, 1, nx); l <= clamp (idx + 1, 1, nx); ++l) {
                for (long m = clamp (jdy - 1, 1, ny); m <= clamp (jdy + 1, 1, ny); ++m) {
                    float cabs = fabs (curv[l][m]);
                    if (cabs > cmax) {
                        cmax = cabs;
                        imax = l;
                        jmax = m;
                    }
                }
            }

            /* Apply criteria to filter out spurious corner sequences
            and add non-spurious corners to the chain */
            if (cmax > -1 && fabs (curv[imax][jmax]) > 0.001 &&
                sgn (current->chain[k - 1].curv) == sgn (curv[imax][jmax])) {
                node_add_to_chain (current, k, imax, jmax, curv[imax][jmax]);
            } else {
                /* Delete chain */
                list_delete (chains, current);
            }
        }
    }
    /* Remove chains that consists only of one position
    FIXME: Is this necessary ? */
    // filter_spurious_chains (&chains, n_iter);

    /* FIXME: Don't know if it really works as intented... */
    least_squares_approx (chains, t, q, nx, ny, v);

    if (DEBUG) {
        /* TODO: visualize the corner angles usign Lemma 1 from the paper */
        write_pgm (v, nx, ny, "/home/danielg/uni/thesis/.tmp/corners.pgm", NULL);
        /* Print chains to image for better visualisation*/
        float **chain_vis;
        alloc_matrix (&chain_vis, nx + 2, ny + 2);

        for (long i = 0; i < nx; ++i) {
            for (long j = 0; j < ny; ++j) {
                chain_vis[i][j] = 0;
            }
        }

        node_ptr current = list_head (chains);
        while (current != NULL) {
            chain_ptr c = node_chain (current);
            for (long k = 0; k < n_iter; ++k) {
                chain_vis[c[k].x][c[k].y] = 255;
            }
            current = node_next (current);
        }

        write_pgm (chain_vis, nx, ny, "/home/danielg/uni/thesis/.tmp/chains.pgm", NULL);
    }

    disalloc_vector (t, n_iter);
    disalloc_matrix (curv, nx + 2, ny + 2);

    return chains;
}


list_ptr test_detection (char *in, char *out) {
    float **u, **f;
    long nx, ny;
    read_pgm_and_allocate_memory (in, &nx, &ny, &u);
    alloc_matrix (&f, nx + 2, ny + 2);
    list_ptr chains = amss_corner_detection (u, nx, ny, 0.1, 1, 20, 1, 0.05f, f);
    if (out != NULL)
        write_pgm (f, nx, ny, out, NULL);
    return chains;
}
