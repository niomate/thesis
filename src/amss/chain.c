#include "chain.h"
#include "../utils.h"
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
    Add new node to beginning of list (is faster than pushing it to the back)
*/
void push_chain

(struct node **head, /* head of the list */
 long x,
 long y,
 float curv,
 long t)

{
    struct node *new = (struct node *)malloc (sizeof (struct node));
    struct corner *chain;

    chain = malloc (t * sizeof (struct corner));

    for (long i = 0; i < t; ++i) {
        chain[i] = (struct corner){ -1, -1, INFINITY };
    }

    chain[0] = (struct corner){ x, y, curv };

    new->chain = chain;
    new->error = 0;
    new->slope = 0;
    new->angle = 0;
    new->x0 = -1;
    new->y0 = -1;
    new->next = *head;
    *head = new;
}


void print_list (struct node **head, long t) {
    struct node *current = *head;

    do {
        print_chain (current->chain, t);
        // printf (", cornerness: %.2f, a:%.2f, b:%.2f\n", current->error, current->a, current->b);
        current = current->next;
    } while (current != NULL);
}


void print_chain (struct corner *chain, long t) {
    for (long i = 0; i < t; ++i) {
        printf ("{index: %ld, x: %ld, y: %ld}->", i, chain[i].x, chain[i].y);
    }
}

void remove_chain (struct node **head, struct node *delete) {
    struct node *current = *head;
    struct node *previous = NULL;
    /* Find node with address */
    while (current != delete) {
        previous = current;
        current = current->next;
    }
    /* Adress is not in list */
    if (current == NULL) {
        return;
    }
    /* Delete head */
    if (previous == NULL) {
        *head = (*head)->next;
        current = *head;
    } else {
        /* Delete arbitrary node */
        previous->next = current->next;
        current = current->next;
    }
}

int cmpfunc (const void *a, const void *b) {
    float x = *(float *)a;
    float y = *(float *)b;
    if (x > y)
        return 1;
    else if (x < y)
        return -1;
    else
        return 0;
}


void cornerness_quantile (struct node **head, long t, long n_corners, float q) {
    float cornerness[n_corners];
    long i = 0;
    float threshold = 0;
    struct node *current = *head;
    struct node *previous = NULL;

    /* TODO: Figure out how to filter out faulty sequences! */
    for (struct node *current = *head; current != NULL; current = current->next) {
        // if (fabs (current->error) == 0) {
        //     continue;
        // } else {
        cornerness[i++] = current->error;
        // }
    }

    /* Sort cornerness in ascending order to compute quantile */
    qsort (cornerness, n_corners, sizeof (float), cmpfunc);
    // for (int i = 0; i < n_corners; ++i) {
    //     printf ("%f, ", cornerness[i]);
    // }
    // printf ("\n");
    long qindex = (long)(q * n_corners);
    // printf ("Quantile corner_index is %ld.\n", qindex);
    threshold = cornerness[qindex];

    // printf ("Threshold is %.2f\n", threshold);
    while (current != NULL) {
        // printf ("Error: %f, Slope: %f, Angle: %f\n", current->error, current->slope, current->angle);
        if (current->error > threshold) {
            /* Error too large, remove node from list */
            if (previous != NULL) {
                previous->next = current->next;
                current = current->next;
            } else { /* current == *head */
                *head = (**head).next;
                current = *head;
            }
        } else {
            previous = current;
            current = current->next;
        }
    }
}
