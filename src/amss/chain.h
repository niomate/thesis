#ifndef chain_h__
#define chain_h__

#include <math.h>

struct pixel;
struct node;
struct list;

typedef struct pixel pixel_t;
typedef struct node node_t;
typedef struct list list_t;
typedef pixel_t *chain_ptr;
typedef node_t *node_ptr;
typedef list_t *list_ptr;
typedef void (*iterator_func) (node_ptr, int);

/* Struct for a single cornertip */
struct pixel {
    long x;
    long y;
    float curv;
};

static const pixel_t EMPTY_PIXEL = (pixel_t){ -1, -1, NAN };

/* Struct for a node in a linked list */
struct node {
    chain_ptr chain;
    pixel_t corner_tip;
    float error;
    float slope;
    float angle;
    node_ptr next;
};

/* Struct for a singly linked list */
struct list {
    node_ptr head;
    int size;
    int chain_length;
};

/* Constructor functions */
list_ptr new_list (int);
node_ptr new_node (int);
chain_ptr new_chain (int);
pixel_t new_pixel (long, long, float);

/* Desctructor functions */
void list_destroy (list_ptr); 

/* Functions for list interaction */
node_ptr list_append_new (list_ptr);
node_ptr list_get (list_ptr, int);
node_ptr list_head (list_ptr);
void list_delete (list_ptr, node_ptr);
int list_size (list_ptr);

/* Print utilities */
void print_list (list_ptr);
void print_node (node_ptr, int);
void print_chain (chain_ptr, int);
void print_pixel (pixel_t);

#endif
