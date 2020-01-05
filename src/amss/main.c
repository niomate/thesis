#include "../utils.h"
#include "amss_corner_detection.h"
#include <getopt.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>


static int DEBUG = 1;

int main (int argc, char *argv[]) {
    // char *filename = "images/binary/angles/angle115-2.pgm";
    // char *output = "output.pgm";
    // test_detection (filename, output);
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