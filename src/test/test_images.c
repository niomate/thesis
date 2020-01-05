#include "../amss/amss_corner_detection.h"
#include "../amss/chain.h"
#include <assert.h>
#include <dirent.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Used in amss_corner_detection to toggle debug messages */
extern int DEBUG = 0;

/* To get M_PI from math.h */
#define _USE_MATH_DEFINES
#define TEST_SUCCESS (1)
#define TEST_FAILURE (0)

/* Expected the detected corner to be around the centre of the image
because of the way the test images were generated */
static const int EXPECTED_X = 256;
static const int EXPECTED_Y = 256;
static const float ANGLE_EPSILON = 0.01;
static const long POS_EPSILON = 3;
static const char *TEST_DIR = "/home/danielg/uni/thesis/images/binary/angles/";

int run_single_test (char *);
void run_all_tests ();

int run_single_test (char *image_name) {
    printf ("Running test for %s...\n", image_name);
    char path[256];
    char outpath[256] = "/home/danielg/uni/thesis/output/";

    /* Extract angle xxx from test image name "anglexxx-y.pgm" */
    char angle_c[3] = { image_name[5], image_name[6], image_name[7] };
    float expected_angle = atof (angle_c) * M_PI / 180;

    strcpy (path, TEST_DIR);
    strcat (path, image_name);
    strcat (outpath, image_name);
    struct node *results = test_detection (path, outpath);

    for (struct node *current = results; current != NULL; current = current->next) {
        if (fabs (current->x0 - EXPECTED_X) <= POS_EPSILON && fabs (current->y0 - EXPECTED_Y) <= POS_EPSILON) {
            /* Check if computed angle matches the actual angle */
            float diff_angle = fabs (current->angle - expected_angle);
            // printf ("Computed angle: %f, expected angle: %f\n", current->angle, expected_angle);
            if (diff_angle < ANGLE_EPSILON)
                return TEST_SUCCESS;
            // else {
            //     printf ("\tDifference in angle: %f\n", diff_angle);
            //     return TEST_FAILURE;
            //     continue;
            // }
        }
    }
    return TEST_FAILURE;
}


void run_all_tests () {
    DIR *p_dir = opendir (TEST_DIR);
    struct dirent *entry;
    int total = 0;
    int failed = 0;

    if (p_dir == NULL) {
        printf ("Cannot open directory %s.\n", TEST_DIR);
        return;
    }

    while ((entry = readdir (p_dir)) != NULL) {
        if (strcmp (entry->d_name, ".") == 0 || strcmp (entry->d_name, "..") == 0)
            continue;
        ++total;
    }

    char **failed_tests = malloc (total * sizeof (char *));
    for (int i = 0; i < total; ++i) {
        failed_tests[i] = malloc (256 * sizeof (char));
    }

    /* Reopen directory since it is mutated by the read operation */
    p_dir = opendir (TEST_DIR);
    while ((entry = readdir (p_dir)) != NULL) {
        if (strcmp (entry->d_name, ".") == 0 || strcmp (entry->d_name, "..") == 0)
            continue;
        if (run_single_test (entry->d_name) == TEST_SUCCESS) {
            printf ("Success!\n");
        } else {
            printf ("Failed!\n");
            failed_tests[failed++] = entry->d_name;
        }
    }

    if (failed == 0) {
        printf ("Passed all tests. Good job!\n");
    } else {
        printf ("Passed %d of %d tests.\n", total - failed, total);
        printf ("Failed %d tests.\n", failed);
        printf ("Following tests did not pass: \n");
        for (int i = 0; i < failed; ++i) {
            printf ("\t%s\n", failed_tests[i]);
        }
    }
    closedir (p_dir);
}

int main () {
    printf ("Running tests for edge detection.\nTest directory: \n  %s\n", TEST_DIR);
    run_all_tests ();
}