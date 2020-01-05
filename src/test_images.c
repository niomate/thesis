#include "amss/amss_corner_detection.h"
#include "amss/chain.h"
#include <assert.h>
#include <dirent.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TEST_SUCCESS (1)
#define TEST_FAILURE (0)

static const int EXPECTED_X = 256;
static const int EXPECTED_Y = 256;

extern int DEBUG = 0;

const char *TEST_DIR = "/home/danielg/uni/thesis/images/binary/angles/";


int run_single_test (char *);
void run_all_tests ();

int run_single_test (char *image_name) {
    printf ("Running test for %s...", image_name);
    char path[256];
    char outpath[256] = "/home/danielg/uni/thesis/output/";
    strcpy (path, TEST_DIR);
    strcat (path, image_name);
    strcat (outpath, image_name);
    struct node *results = test_detection (path, outpath);

    for (struct node *current = results; current != NULL; current = current->next) {
        // printf ("%ld, %ld\n", current->x0, current->y0);
        if (fabs (current->x0 - EXPECTED_X) <= 3 && fabs (current->y0 - EXPECTED_Y) <= 3) {
            return TEST_SUCCESS;
        }
    }
    return TEST_FAILURE;
}


void run_all_tests (const char *dir) {
    DIR *p_dir = opendir (dir);
    struct dirent *entry;

    if (p_dir == NULL) {
        printf ("Cannot open directory %s.\n", dir);
        return;
    }

    int total = 0;
    int success = 0;

    while ((entry = readdir (p_dir)) != NULL) {
        if (strcmp (entry->d_name, ".") == 0 || strcmp (entry->d_name, "..") == 0)
            continue;
        if (run_single_test (entry->d_name) == TEST_SUCCESS) {
            printf ("Success!\n");
            ++success;
        } else {
            printf ("Failed!\n");
        }
        ++total;
    }

    if (success == total) {
        printf ("Passed all tests. Good job!\n");
    } else {
        printf ("Passed %d of %d tests.\n", success, total);
        printf ("Failed %d tests.\n", total - success);
    }

    closedir (p_dir);
}

int main () {
    printf ("Running tests for edge detection.\nTest directory: \n  %s\n", TEST_DIR);
    run_all_tests (TEST_DIR);
}