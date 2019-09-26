#ifndef utils_h__
#define utils_h__

float quantile (float* u, float q, long n);

float quantile_2D (float **u, float q, long nx, long ny);

float MSE (float **u, float **v, long nx, long ny);

float MSE_filenames (char *im1, char *im2);

void alloc_vector (float **vector, long n);

void alloc_matrix (float ***matrix, long nx, long ny);

void alloc_cubix (float ****cubix, long nx, long ny, long nz);

void disalloc_vector (float *vector, long n);

void disalloc_matrix (float **matrix, long nx, long ny);

void disalloc_cubix (float ***cubix, long nx, long ny, long nz);

#endif // utils_h__