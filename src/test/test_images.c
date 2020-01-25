#include "../amss/amss_corner_detection.h"
#include "../amss/chain.h"
#include <assert.h>
#include <dirent.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* To get M_PI from math.h */
#define _USE_MATH_DEFINES

/* Expected the detected corner to be around the centre of the image
because of the way the test images were generated */
static const int EXPECTED_X = 64;
static const int EXPECTED_Y = 64;
static const float ANGLE_EPSILON = 0.01;
static const long POS_EPSILON = 2;
static const float ERROR_EPSILON = 0.01;
static const char *TEST_DIR = "/home/danielg/uni/thesis/images/binary/angles/";
static const char *TEST_LOG = "/home/danielg/uni/thesis/test.log";

int run_single_test_case (char *);
void run_all_tests ();

int test_position (list_ptr results) {
    FILE *log = fopen (TEST_LOG, "a");
    fprintf (log, "\t[test_position]\n");
    int success = 0;
    for (node_ptr current = list_head (results); current != NULL; current = current->next) {
        long x = current->corner_tip.x;
        long y = current->corner_tip.y;
        float pos_error = 0.5 * (powf (x - EXPECTED_X, 2.0f) + powf (y - EXPECTED_Y, 2.0f));

        if (pos_error < POS_EPSILON) {
            success = 1;
            break;
        }
    }
    if (!success)
        fprintf (log, "\tTest failed!\n");
    else
        fprintf (log, "\tTest passed!\n");
    fclose (log);
    return success;
}

int test_angle (list_ptr results, float expected_angle) {
    FILE *log = fopen (TEST_LOG, "a");
    fprintf (log, "\t[test_angle]\n");
    int success = 0;
    for (node_ptr current = list_head (results); current != NULL; current = current->next) {
        long x = current->corner_tip.x;
        long y = current->corner_tip.y;
        float pos_error = 0.5 * (powf (x - EXPECTED_X, 2.0f) + powf (y - EXPECTED_Y, 2.0f));
        if (pos_error < POS_EPSILON) {
            /* Check if computed angle matches the actual angle */
            float current_angle = current->angle;
            float diff_angle = fabs (current_angle - expected_angle);
            fprintf (log, "\t\tComputed angle: %f, expected angle: %f\n", current_angle, expected_angle);
            if (diff_angle < ANGLE_EPSILON) {
                success = 1;
                break;
            } else {
                fprintf (log, "\t\tDifference in angle: %f\n", diff_angle);
                fprintf (log, "\t\tTrying next corner\n");
                continue;
            }
        }
    }
    if (!success)
        fprintf (log, "\tTest failed!\n");
    else
        fprintf (log, "\tTest passed!\n");
    fclose (log);
    return success;
}

int test_error (list_ptr results) {
    FILE *log = fopen (TEST_LOG, "a");
    fprintf (log, "\t[test_residual_error]\n");
    int success = 0;
    for (node_ptr current = list_head (results); current != NULL; current = current->next) {
        long x = current->corner_tip.x;
        long y = current->corner_tip.y;
        float pos_error = 0.5 * (powf (x - EXPECTED_X, 2.0f) + powf (y - EXPECTED_Y, 2.0f));
        if (pos_error < POS_EPSILON) {
            if (current->error < ERROR_EPSILON) {
                success = 1;
                break;
            } else {
                fprintf (log, "\t\tResidual error: %f\n", current->error);
                fprintf (log, "\t\tTrying next corner\n");
                continue;
            }
        }
    }
    if (!success)
        fprintf (log, "\tTest failed!\n");
    else
        fprintf (log, "\tTest passed!\n");
    fclose (log);
    return success;
}

int run_single_test_case (char *image_name) {
    FILE *log = fopen (TEST_LOG, "a");
    fprintf (log, "[Running tests for %s...]\n", image_name);
    fclose (log);
    char path[256];
    char outpath[256] = "/home/danielg/uni/thesis/output/";

    /* Extract angle xxx from test image name "anglexxx-y.pgm" */
    char angle_c[3] = { image_name[5], image_name[6], image_name[7] };
    float expected_angle = atof (angle_c);
    float expected_slope = 1.0f / sqrt (tan (expected_angle * (M_PI / 180) * 0.5));
    printf ("Expecting a slope of %f\n", expected_slope);

    /* TODO: Insecure! */
    strcpy (path, TEST_DIR);
    strcat (path, image_name);
    strcat (outpath, image_name);

    list_ptr results = test_detection (path, outpath);

    int position_result = test_position (results);
    int angle_result = test_angle (results, expected_angle);
    int error_result = test_error (results);

    return position_result + angle_result + error_result;
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

    char **failed_tests = malloc (total * sizeof (*failed_tests));
    for (int i = 0; i < total; ++i) {
        failed_tests[i] = malloc (256 * sizeof (char));
    }

    /* Reopen directory since it is mutated by the read operation */
    p_dir = opendir (TEST_DIR);
    while ((entry = readdir (p_dir)) != NULL) {
        if (strcmp (entry->d_name, ".") == 0 || strcmp (entry->d_name, "..") == 0)
            continue;

        //if (strcmp (entry->d_name, "angle090-0.pgm") != 0) {
            //continue;
        //}
        printf ("[%s]: ", entry->d_name);
        int n = 0;
        if ((n = run_single_test_case (entry->d_name)) == 3) {
            printf ("Success!\n");
        } else {
            printf ("Failed %d/3 tests!\n", 3 - n);
            failed_tests[failed++] = entry->d_name;
        }
    }

    FILE *log = fopen (TEST_LOG, "a");

    if (failed == 0) {
        fprintf (log, "Passed all tests. Good job!\n");
    } else {
        fprintf (log, "Passed %d of %d tests.\n", total - failed, total);
        fprintf (log, "Failed %d tests.\n", failed);
        fprintf (log, "Following tests did not pass: \n");
        for (int i = 0; i < failed; ++i) {
            fprintf (log, "\t%s\n", failed_tests[i]);
        }
    }

    fclose (log);
    closedir (p_dir);
}

int main () {
    FILE *log = fopen (TEST_LOG, "w");
    fprintf (log, "Running tests for edge detection.\nTest directory: \n  %s\n", TEST_DIR);
    printf ("Running tests for edge detection.\nTest directory: \n  %s\n", TEST_DIR);
    fclose (log);
    run_all_tests ();
}
