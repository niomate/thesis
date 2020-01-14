#ifndef chain_h__
#define chain_h__

#include <math.h>
#include <stdbool.h>

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

/* Functions for list interaction */
list_ptr list_new (int);
node_ptr list_insert (list_ptr, node_ptr);
node_ptr list_get (list_ptr, int);
node_ptr list_head (list_ptr);
bool list_end (node_ptr);
void list_delete (list_ptr, node_ptr);
void list_foreach (list_ptr, iterator_func);
int list_size (list_ptr);
int list_chain_length (list_ptr);

/* Functions for node interaction */
node_ptr node_new (int);
node_ptr node_next (node_ptr);
chain_ptr node_chain (node_ptr);
pixel_t node_chain_get_pixel (node_ptr, long);
void node_add_to_chain (node_ptr, long, long, long, float);

/* Functions for chain type */
chain_ptr chain_new (int);

/* Functions for pixel type */
pixel_t pixel_new (long, long, float);

/* Print utilities */
void print_list (list_ptr);
void print_node (node_ptr, int);
void print_chain (chain_ptr, int);
void print_pixel (pixel_t);

void list_destroy (list_ptr); 
#endif
