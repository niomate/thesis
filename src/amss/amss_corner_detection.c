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

static int DEBUG = false;

/*
    Compute the differential operator L(u) as stated in the paper ov Alvarez
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
            u[i][j] = (u[i][j] - min) * 255.0 / (max - min);
        }
    }
}

void corner_estimation (struct node **chains, float *t, float q, long n_iter, long nx, long ny, float **corners) {
    /* Least squares approximation */
    float **M; /* System matrix */
    float *b;  /* Value vector */
    float *w;  /* Weight vector, basically just an array of ones */
    float *x;  /* Solution of linear system */

    /* !!! CARE: Iteration in least squares starts with 1 !!! */
    alloc_matrix (&M, n_iter + 1, 3);
    alloc_vector (&b, n_iter + 1);
    alloc_vector (&w, n_iter + 1);
    alloc_vector (&x, 3);

    for (long k = 0; k < n_iter + 1; ++k) {
        /* Equal weights for least squares computation */
        w[k] = 1;

        M[k][1] = 1;
        M[k][2] = 1;
    }

    for (long k = 0; k < n_iter; ++k) {
        M[k + 1][1] = powf (t[k], 0.75f) - powf (t[0], 0.75f);
    }

    long n_corners_before = 0;
    long n_corners_after = 0;

    for (struct node *current = *chains; current != NULL; current = current->next) {
        struct corner *chain = current->chain;
        /* Initialise value vector */
        for (long k = 0; k < n_iter; ++k) {
            /* Distance to first point of corner chain */
            b[k + 1] = sqrt (powf (chain[k].x - chain[0].x, 2.0) + powf (chain[k].y - chain[0].y, 2.0));
        }

        /* Approximate line using least squares approximation */
        solve_normal_equations (n_iter, 2, M, b, w, x);

        /* Need to correct slope from estimation as stated in the paper */
        // current->slope = 0.071917f + 1.029484f * x[1];
        current->slope = x[1];
        current->angle = 2 * atan2f (1, powf (current->slope, 2.0f));

        /* Calculate cornerness a.k.a approximation error */
        float error = 0;
        float acc;

        for (long k = 1; k < n_iter + 1; ++k) {
            acc = b[k] - (current->slope * M[k][1] + x[2]);
            error += pow (acc, 2.0f);
        }

        error /= n_iter;

        printf ("Error with a=%f, b=%f: %f\n", x[1], x[2], error);
        current->error = error;

        /* TODO: make it prettier */
        /* TODO: Filter out chains with slope smaller than epsilon ? */
        float t0x, t0y, tmaxx, tmaxy;

        t0x = current->chain[0].x;
        t0y = current->chain[0].y;
        tmaxx = current->chain[n_iter - 1].x;
        tmaxy = current->chain[n_iter - 1].y;

        /* Reconstruct initial corner position by extrapolation using
        the corner bisector unit vector and slope */

        float norm_inv = 1 / (sqrt (powf (tmaxx - t0x, 2.0f) + powf (tmaxy - t0y, 2.0f)) + FLT_EPSILON);
        float helper = current->slope * powf (1.33333333f * t[0], 0.75f) * norm_inv;

        long x0 = (long)(t0x - helper * (tmaxx - t0x));
        long y0 = (long)(t0y - helper * (tmaxy - t0y));

        current->x0 = x0;
        current->y0 = y0;

        ++n_corners_before;
    }

    /* Only keep quantile of corner chains */
    printf ("Found %ld corner sequences.\n", n_corners_before);
    printf ("Applying quantile to corner chains...\n");

    /* TODO: Figure out quantile parameter */
    // cornerness_quantile (chains, n_iter, n_corners_before, q);

    for (long i = 1; i <= nx; ++i) {
        for (long j = 1; j <= ny; ++j) {
            corners[i][j] = 0;
        }
    }

    for (struct node *current = *chains; current != NULL; current = current->next) {
        ++n_corners_after;
        corners[current->x0][current->y0] = 255.0;
    }

    printf ("Corners left after quantile application: %ld .\n", n_corners_after);

    disalloc_matrix (M, n_iter + 1, 3);
    disalloc_vector (b, n_iter + 1);
    disalloc_vector (w, n_iter + 1);
    disalloc_vector (x, 3);
}

/*
    Detect corners in the given image using the AMSS.
*/
void amss_corner_detection

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
    long i, j, k, l, m;                  /* loop variables */
    long n_iter = t_max / step + 1;      /* number of total iterations */

    /* variables for corner tracking */
    long imax, jmax;
    float cmax;
    float dx, dy, norm_inv;
    long idx, jdy;
    int n_corners = 0;

    /* Linked list to store corner chains in */
    struct node *chains = NULL;

    if (DEBUG)
        printf ("Number of iterations: %ld\n", n_iter);

    alloc_vector (&t, n_iter);
    alloc_matrix (&curv, nx + 2, ny + 2);

    printf ("Initial scale: %ld", t_0);

    t[0] = t_0;

    for (float t = 0; t < t_0; t += ht) {
        amss (ht, nx, ny, h, h, u);
    }

    curvature (u, nx, ny, h, curv);
    dummies (curv, nx, ny);

    /* Find initial set of corners */
    for (i = 1; i <= nx; i++) {
        for (j = 1; j <= ny; j++) {
            /* if u[i][j] is bigger (or smaller) than all of its 8 neighbours
            then it is considered a local extremum */
            if ((fabs (curv[i][j]) > fabs (curv[i + 1][j])) && (fabs (curv[i][j]) > fabs (curv[i - 1][j])) &&
                (fabs (curv[i][j]) > fabs (curv[i][j + 1])) && (fabs (curv[i][j]) > fabs (curv[i][j - 1])) &&
                (curv[i + 1][j] * curv[i - 1][j] * curv[i][j + 1] * curv[i][j - 1] != 0.0)) {
                ++n_corners;
                push_chain (&chains, i, j, curv[i][j], n_iter);
            }
        }
    }

    printf ("Found %d extrema at the initial scale\n", n_corners);

    /* Track corners across time scales */
    for (k = 1; k < n_iter; ++k) {
        printf ("Iteration %ld\n", k);
        n_corners = 0;
        t[k] = t[k - 1] + step;

        /* First: calculate amss at time level t[k] in-place */
        for (float t = 0; t < step; t += ht) {
            amss (ht, nx, ny, h, h, u);
        }
        curvature (u, nx, ny, h, curv);
        dummies (curv, nx, ny);

        /* Find new maximum for each chain at new timescale */
        for (struct node *current = chains; current != NULL; current = current->next) {
            n_corners++;

            /* Position of corner at last scale */
            i = current->chain[k - 1].x;
            j = current->chain[k - 1].y;

            assert (i >= 1);
            assert (j >= 1);

            /* Compute first order derivatives by means of central differences */
            dx = (u[i - 1][j] - u[i + 1][j]) * two_h_inv;
            dy = (u[i][j - 1] - u[i][j + 1]) * two_h_inv;

            /* Norm gradient to 1 */
            norm_inv = 1 / (sqrt (dx * dx + dy * dy) + FLT_EPSILON);
            dx *= norm_inv;
            dy *= norm_inv;

            /* New locations for extrema search */
            idx = (long)round (i + dx);
            jdy = (long)round (j + dy);

            /* Mirroring boundary conditions (correct ?) */
            if (idx < 0) {
                idx = 1 - idx;
            } else if (idx > nx) {
                idx = 2 * nx + 1 - idx;
            }

            if (jdy < 0) {
                jdy = 1 - jdy;
            } else if (jdy > ny) {
                jdy = 2 * ny + 1 - jdy;
            }

            imax = jmax = -1;
            cmax = -1;

            /* Search for new curvature extremum in 8-neighbourhood around pixel (idx, jdy) */
            for (l = idx - 1; l <= idx + 1; ++l) {
                for (m = jdy - 1; m <= jdy + 1; ++m) {
                    if (abs (curv[l][m]) > cmax) {
                        cmax = abs (curv[l][m]);
                        imax = l;
                        jmax = m;
                    }
                }
            }

            /* Apply criteria to filter out spurious corner sequences
            and add non-spurious corners to the chain */
            if (cmax > -1 && fabs (curv[imax][jmax]) > 0.001 &&
                sgn (current->chain[k - 1].curv) == sgn (curv[imax][jmax])) {
                /* Add to chain */
                current->chain[k] = (struct corner){ imax, jmax, curv[imax][jmax] };
            } else {
                /* Delete chain */
                remove_chain (&chains, current);
            }
        }
        printf ("Found %d corners in iteration %ld.\n", n_corners, k);
        char name[256];
        sprintf (name, "/home/danielg/uni/thesis/.tmp/curv%ld.pgm", k);
        normalize_range (curv, nx, ny);
        write_pgm (curv, nx, ny, name, NULL);
    }

    corner_estimation (&chains, t, q, n_iter, nx, ny, v);

    write_pgm (v, nx, ny, "/home/danielg/uni/thesis/.tmp/corners.pgm", NULL);

    if (DEBUG) {
        /* Print chains to image for better visualisation*/
        float **chain_vis;
        alloc_matrix (&chain_vis, nx + 2, ny + 2);

        struct node *current = chains;

        for (i = 0; i < nx; ++i) {
            for (j = 0; j < ny; ++j) {
                chain_vis[i][j] = 0;
            }
        }

        while (current != NULL) {
            struct corner *c = current->chain;
            // for (k = 0; k < n_iter; ++k) {
            //     chain_vis[c[k].x][c[k].y] = 255;
            // }
            chain_vis[c[0].x][c[0].x] = 255.0;
            current = current->next;
        }

        write_pgm (chain_vis, nx, ny, "/home/danielg/uni/thesis/.tmp/chains.pgm", NULL);
    }

    disalloc_vector (t, n_iter);
    disalloc_matrix (curv, nx + 2, ny + 2);
}


int main (int argc, char *argv[]) {
    char *filename; /* filename for image; gets read from command line */
    float **f;      /* original image */
    float **mask;   /* mask: locations of corners */
    long nx, ny;    /* dimensions of the image */

    float ht;           /* timestep size */
    long t_0;           /* start scale for AMSS */
    long t_max, t_step; /* initial scale and scale step size */

    float quant; /* quantile parameter for extrema selection */

    /* variables for command line parsing */
    int option_index = 0; /* contains corner_index of current option */
    int opt;              /* contains current option */

    static struct option long_options[] = { { "debug", no_argument, &DEBUG, true },
                                            { "ht", required_argument, 0, 't' },
                                            { "t_max", required_argument, 0, 'm' },
                                            { "t_step", required_argument, 0, 's' },
                                            { "quantile", required_argument, 0, 'q' },
                                            { "t_0", required_argument, 0, 'i' },
                                            { 0, 0, 0, 0 } };


    /* default options */
    ht = 0.1;
    t_0 = 0;
    t_max = 1000;
    t_step = 100;
    quant = 0.05;

    /* read command line arguments */
    while ((opt = getopt_long (argc, argv, "t:m:s:q:i:", long_options, &option_index)) != -1) {
        switch (opt) {
        case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0)
                break;
            printf ("option %s", long_options[option_index].name);
            if (optarg)
                printf (" with arg %s", optarg);
            printf ("\n");
            break;
        case 'd':
            DEBUG = true;
            break;
        case 't':
            ht = atof (optarg);
            break;
        case 'm':
            t_max = atol (optarg);
            break;
        case 's':
            t_step = atol (optarg);
            break;
        case 'q':
            quant = atof (optarg);
            break;
        case 'i':
            t_0 = atol (optarg);
            break;
        default:
            fprintf (stderr,
                     "Usage: %s [-t <timestep size>] "
                     "[-m <max scale>] [-s <scalestep size>]"
                     "[-q <quantile>] [-i <initial scale>] [file]\n",
                     argv[0]);
            exit (EXIT_FAILURE);
        }
    }

    if (optind > argc - 1) {
        fprintf (stderr, "No input image specified. Exiting...\n");
        return -1;
    }

    filename = argv[optind];


    read_pgm_and_allocate_memory (filename, &nx, &ny, &f);
    alloc_matrix (&mask, nx + 2, ny + 2);
    amss_corner_detection (f, nx, ny, ht, t_0, t_max, t_step, quant, mask);

    return 1;
}