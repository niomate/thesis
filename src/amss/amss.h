#ifndef amss_h__
#define amss_h__

/*
 affine morphological scale-space, explicit scheme
*/
void amss (float ht, long nx, long ny, float hx, float hy, float **u);

void read_string (char *v);

void read_long (long *v);

void read_float (float *v);

void read_pgm_and_allocate_memory (const char *file_name, long *nx, long *ny, float ***u);

void comment_line (char *comment, char *lineformat, ...);

void write_pgm (float **u, long nx, long ny, char *file_name, char *comments);

void dummies (float **u, long nx, long ny);

void analyse_grey (float **u, long nx, long ny, float *min, float *max, float *mean, float *std);

float sgn(float x);

#endif // amss_h__