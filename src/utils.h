#ifndef utils_h__
#define utils_h__

void alloc_vector (float **, long);

void alloc_vector_long (long **, long);

void alloc_matrix (float ***, long, long);

void alloc_cubix (float ****, long, long, long);

void disalloc_vector (float *, long);

void disalloc_matrix (float **, long, long);

void disalloc_cubix (float ***, long, long, long);

void read_string (char *);

void read_long (long *);

void read_float (float *);

void read_pgm_and_allocate_memory (const char *, long *, long *, float ***);

void comment_line (char *, char *, ...);

void write_pgm (float **, long, long, char *, char *);

void dummies (float **, long, long);

float quantile (float *, long, float);

float quantile_2D (float **, long, long, float);

float sgn (float);

void analyse_grey (float **, long, long, float *, float *, float *, float *);

void mask (float **, float **, long, long, int);

long clamp (long, long, long);

#endif // utils_h__
