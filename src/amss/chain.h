#ifndef chain_h__
#define chain_h__

/* Struct for a single cornertip */
struct corner {
    long x;
    long y;
    float curv;
};

/* Struct for a node in a linked list */
struct node {
    struct corner *chain;
    float error;
    float slope;
    float angle;
    long x0;
    long y0;
    struct node *next;
};

void push_chain (struct node **head, long x, long y, float curv, long t);

void print_list (struct node **head, long t);

void print_chain (struct corner *chain, long t);

void print_complete (struct node **head, long t);

void remove_chain (struct node **head, struct node *delete);

void cornerness_quantile (struct node **head, long t, long n_corners, float q);


#endif