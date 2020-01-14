#include "chain.h"
#include "../utils.h"
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>


list_ptr list_new (int chain_length) {
    list_ptr list = (list_ptr)malloc (sizeof (list_t));
    list->head = NULL;
    list->size = 0;
    list->chain_length = chain_length;
    return list;
}


node_ptr list_insert (list_ptr list, node_ptr node) {
    if (!node) {
        node = node_new (list->chain_length);
    }

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


void list_foreach (list_ptr list, iterator_func func) {
    for (node_ptr current = list->head; current != NULL; current = current->next) {
        func (current, list_size (list));
    }
}


int list_size (list_ptr list) {
    return list->size;
}


int list_chain_length (list_ptr list) {
    return list->chain_length;
}


node_ptr node_new (int chain_length) {
    node_ptr n = (node_ptr)malloc (sizeof (node_t));
    n->error = NAN;
    n->slope = NAN;
    n->angle = NAN;
    n->next = NULL;
    n->corner_tip = EMPTY_PIXEL;
    n->chain = chain_new (chain_length);
    return n;
}


node_ptr node_next (node_ptr node) {
    return node->next;
}


chain_ptr node_chain (node_ptr node) {
    return node->chain;
}


void node_add_to_chain (node_ptr node, long k, long x, long y, float curv) {
    node->chain[k] = pixel_new (x, y, curv);
}


chain_ptr chain_new (int size) {
    chain_ptr c = malloc (size * sizeof (pixel_t));
    for (int i = 0; i < size; ++i) {
        c[i] = EMPTY_PIXEL;
    }
    return c;
}


pixel_t pixel_new (long x, long y, float curv) {
    return (pixel_t){ .x = x, .y = y, .curv = curv };
}


void print_list (list_ptr list) {
    printf ("List[");
    list_foreach (list, print_node);
    printf ("]\n");
}


void print_node (node_ptr node, int chain_length) {
    printf ("Node [");
    print_chain (node->chain, chain_length);
    print_pixel (node->corner_tip);
    printf ("error: %f, slope: %f, angle: %f]->", node->error, node->slope, node->angle);
    printf ("]->");
}


void print_chain (chain_ptr chain, int chain_length) {
    printf ("Chain [");
    for (int i = 0; i < chain_length; ++i) {
        print_pixel (chain[i]);
    }
    printf ("] ");
}

void print_pixel (pixel_t pixel) {
    printf ("Pixel [%ld, %ld, %f] ", pixel.x, pixel.y, pixel.curv);
}

void list_destroy (list_ptr list) {
    node_ptr current;
    while ((current = list->head) != NULL) {
        list->head = current->next;

        free (current->chain);
        free (current);
    }
}
