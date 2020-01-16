#include "chain.h"
#include "../utils.h"
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

/****************************** Constructor functions **********************************/

list_ptr new_list (int chain_length) {
    list_ptr list = (list_ptr)malloc (sizeof (list_t));
    list->head = NULL;
    list->size = 0;
    list->chain_length = chain_length;
    return list;
}


node_ptr new_node (int chain_length) {
    node_ptr n = (node_ptr)malloc (sizeof (node_t));
    n->error = NAN;
    n->slope = NAN;
    n->angle = NAN;
    n->next = NULL;
    n->corner_tip = EMPTY_PIXEL;
    n->chain = new_chain (chain_length);
    return n;
}


chain_ptr new_chain (int size) {
    chain_ptr c = malloc (size * sizeof (pixel_t));
    for (int i = 0; i < size; ++i) {
        c[i] = EMPTY_PIXEL;
    }
    return c;
}


pixel_t new_pixel (long x, long y, float curv) {
    return (pixel_t){ .x = x, .y = y, .curv = curv };
}

/****************************** List desctructor **********************************/

void list_destroy (list_ptr list) {
    node_ptr current;
    while ((current = list->head) != NULL) {
        list->head = current->next;

        free (current->chain);
        free (current);
    }
}

/****************************** List functions **********************************/

node_ptr list_append_new (list_ptr list) {
    node_ptr node = new_node (list->chain_length);
    node->next = list->head;
    list->head = node;
    list->size++;
    return node;
}


node_ptr list_get (list_ptr list, int idx) {
    assert (idx > 0);
    assert (idx < list_size (list));
    node_ptr current = list->head;

    int i = 0;
    while (i < idx) {
        current = current->next;
    }

    return current;
}


node_ptr list_head (list_ptr list) {
    return list->head;
}


bool list_end (node_ptr node) {
    return node->next == NULL;
}


void list_delete (list_ptr list, node_ptr node) {
    node_ptr current = list->head;
    node_ptr prev = NULL;

    while (current != node) {
        prev = current;
        current = current->next;
    }

    if (current == NULL) {
        return;
    }

    if (prev == NULL) {
        /* current == list->head */
        list->head = current->next;
    } else {
        prev->next = current->next;
    }

    list->size--;

    free (current->chain);
    free (current);
}


int list_size (list_ptr list) {
    return list->size;
}

/****************************** Print functions **********************************/

void print_list (list_ptr list) {
    printf ("List[");
    printf ("chain_length: %d, \n", list->chain_length);
    printf ("size: %d, \n", list->size);
    for (node_ptr current = list_head(list); current != NULL; current=current->next) {
      print_node(current, list->chain_length);
    }
    printf ("]\n");
}


void print_node (node_ptr node, int chain_length) {
    printf ("Node [");
    print_chain (node->chain, chain_length);
    print_pixel (node->corner_tip);
    printf ("error: %f, slope: %f, angle: %f]->", node->error, node->slope, node->angle);
    printf ("]\n");
}


void print_chain (chain_ptr chain, int chain_length) {
    printf ("Chain [");
    for (int i = 0; i < chain_length; ++i) {
        print_pixel (chain[i]);
    }
    printf ("]\n");
}

void print_pixel (pixel_t pixel) {
    printf ("Pixel [%ld, %ld] ", pixel.x, pixel.y);
}
