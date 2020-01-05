#ifndef utils_h__
#define utils_h__

void alloc_vector (float **vector, long n);

void alloc_vector_long (long **vector, long n);

void alloc_matrix (float ***matrix, long nx, long ny);

void alloc_cubix (float ****cubix, long nx, long ny, long nz);

void disalloc_vector (float *vector, long n);

void disalloc_matrix (float **matrix, long nx, long ny);

void disalloc_cubix (float ***cubix, long nx, long ny, long nz);

void read_string (char *v);

void read_long (long *v);

void read_float (float *v);

void read_pgm_and_allocate_memory (const char *file_name, long *nx, long *ny, float ***u);

void comment_line (char *comment, char *lineformat, ...);

void write_pgm (float **u, long nx, long ny, char *file_name, char *comments);

void dummies (float **u, long nx, long ny);

float quantile (float *arr, long n, float q);

float quantile_2D (float **u, long nx, long ny, float q);

float sgn (float x);

void analyse_grey (float **u, long nx, long ny, float *min, float *max, float *mean, float *std);

void mask (float **u, float **m, long nx, long ny, int r);

long clamp (long x, long lo, long hi);

#endif // utils_h__