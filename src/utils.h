#ifndef utils_h__
#define utils_h__

#include <stdlib.h>

#ifndef MAX
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#endif

/* Matrix alloc and disalloc utils */
void alloc_vector(float **u, long n);
void alloc_vector_long(long **u, long n);
void alloc_matrix(float ***u, long nx, long ny);
void alloc_cubix(float ****u, long nx, long ny, long nz);
void disalloc_vector(float *u, long n);
void disalloc_matrix(float **u, long nx, long ny);
void disalloc_cubix(float ***u, long nx, long ny, long nz);

/* Read utilities */
void read_string(char *v);
void read_long(long *v);
void read_float(float *v);

/* Image utilities */
void read_pgm_and_allocate_memory(const char *file_name, long *nx, long *ny,
                                  float ***u);
void comment_line(char *comment, char *format, ...);
void write_pgm(float **u, long nx, long ny, char *file_name, char *comment);
void dummies(float **u, long nx, long ny);
void imgcpy(float **src, float **dst, size_t nx, size_t ny);

/* Math utilities */
float sgn(float x);
void analyse_grey(float **u, long nx, long ny, float *min, float *max,
                  float *mean, float *std);
long clamp(long n, long lo, long hi);
int in_circle(long x, long y, long cx, long cy, float radius);
int out_of_bounds(long i, long n);
int in_image(long i, long j, long nx, long ny);
int gauss_circle(float radius);

#endif // utils_h__
