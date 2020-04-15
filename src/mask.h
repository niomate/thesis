#ifndef mask_h__
#define mask_h__

void mask(float **u, long nx, long ny, int radius, float **v);
void sample_mask (float **u, long nx, long ny, int radius, float perc, float **v);

#endif  // mask_h__
