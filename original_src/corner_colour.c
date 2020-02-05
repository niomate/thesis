#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define pi 3.1415927

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*               CORNER DETECTION WITH THE STRUCTURE TENSOR                 */
/*                                                                          */
/*                  (Copyright Joachim Weickert, 12/2005)                   */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/* 
 features:
 - for greyscale and colour images
 - presmoothing at noise scale:  convolution-based, Neumann b.c.
 - presmoothing at integration scale: convolution-based, Dirichlet b.c.
*/


/*--------------------------------------------------------------------------*/

void alloc_vector

     (float **vector,   /* vector */
      long  n)          /* size */

     /* allocates storage for a vector of size n */


{
*vector = (float *) malloc (n * sizeof(float));
if (*vector == NULL)
   {
   printf("alloc_vector: not enough storage available\n");
   exit(1);
   }
return;
}

/*--------------------------------------------------------------------------*/

void alloc_matrix

     (float ***matrix,  /* matrix */
      long  nx,         /* size in x direction */
      long  ny)         /* size in y direction */

     /* allocates storage for matrix of size nx * ny */


{
long i;

*matrix = (float **) malloc (nx * sizeof(float *));
if (*matrix == NULL)
   {
   printf("alloc_matrix: not enough storage available\n");
   exit(1);
   }
for (i=0; i<nx; i++)
    {
    (*matrix)[i] = (float *) malloc (ny * sizeof(float));
    if ((*matrix)[i] == NULL)
       {
       printf("alloc_matrix: not enough storage available\n");
       exit(1);
       }
    }
return;
}

/*--------------------------------------------------------------------------*/

void alloc_cubix

     (float ****cubix,  /* cubix */
      long  nx,         /* size in x direction */
      long  ny,         /* size in y direction */
      long  nz)         /* size in z direction */

     /* allocates storage for cubix of size nx * ny * nz */


{
long i, j;

*cubix = (float ***) malloc (nx * sizeof(float **));
if (*cubix == NULL)
   {
   printf("alloc_cubix: not enough storage available\n");
   exit(1);
   }
for (i=0; i<nx; i++)
    {
    (*cubix)[i] = (float **) malloc (ny * sizeof(float *));
    if ((*cubix)[i] == NULL)
       {
       printf("alloc_cubix: not enough storage available\n");
       exit(1);
       }
    for (j=0; j<ny; j++)
        {
        (*cubix)[i][j] = (float *) malloc (nz * sizeof(float));
        if ((*cubix)[i][j] == NULL)
           {
           printf("alloc_cubix: not enough storage available\n");
           exit(1);
           }
        }
    }
return;
}

/*--------------------------------------------------------------------------*/

void disalloc_vector

     (float *vector,    /* vector */
      long  n)          /* size */

     /* disallocates storage for a vector of size n */

{
free(vector);
return;
}

/*--------------------------------------------------------------------------*/

void disalloc_matrix

     (float **matrix,   /* matrix */
      long  nx,         /* size in x direction */
      long  ny)         /* size in y direction */

     /* disallocates storage for matrix of size nx * ny */

{
long i;
for (i=0; i<nx; i++)
    free(matrix[i]);
free(matrix);
return;
}

/*--------------------------------------------------------------------------*/

void disalloc_cubix

     (float ***cubix,   /* cubix */
      long  nx,         /* size in x direction */
      long  ny,         /* size in y direction */
      long  nz)         /* size in z direction */

     /* disallocates storage for cubix of size nx * ny * nz */

{
long i, j;
for (i=0; i<nx; i++)
 for (j=0; j<ny; j++)
     free(cubix[i][j]);
for (i=0; i<nx; i++)
    free(cubix[i]);
free(cubix);
return;
}

/* ----------------------------------------------------------------------- */

void dummies

     (float **v,        /* image matrix */
      long  nx,         /* size in x direction */
      long  ny)         /* size in y direction */

/* creates dummy boundaries by mirroring */

{
long i, j;  /* loop variables */

for (i=1; i<=nx; i++)
    {
    v[i][0]    = v[i][1];
    v[i][ny+1] = v[i][ny];
    }

for (j=0; j<=ny+1; j++)
    {
    v[0][j]    = v[1][j];
    v[nx+1][j] = v[nx][j];
    }
return;
}

/* ---------------------------------------------------------------------- */

void gauss_conv 

     (float    sigma,     /* standard deviation of Gaussian */
      long     nx,        /* image dimension in x direction */ 
      long     ny,        /* image dimension in y direction */ 
      float    hx,        /* pixel size in x direction */
      float    hy,        /* pixel size in y direction */
      float    precision, /* cutoff at precision * sigma */
      long     bc,        /* type of boundary condition */
                          /* 0=Dirichlet, 1=reflecing, 2=periodic */
      float    **f)       /* input: original image ;  output: smoothed */


/* 
 Gaussian convolution. 
*/


{
long    i, j, p;              /* loop variables */
long    length;               /* convolution vector: 0..length */
float   sum;                  /* for summing up */
float   *conv;                /* convolution vector */
float   *help;                /* row or column with dummy boundaries */
      

/* ------------------------ diffusion in x direction -------------------- */

/* calculate length of convolution vector */
length = (long)(precision * sigma / hx) + 1;
if ((bc != 0) && (length > nx))
   {
   printf("gauss_conv: sigma too large \n"); 
   exit(0);
   }

/* allocate storage for convolution vector */
alloc_vector (&conv, length+1);

/* calculate entries of convolution vector */
for (i=0; i<=length; i++)
    conv[i] = 1 / (sigma * sqrt(2.0 * 3.1415926)) 
              * exp (- (i * i * hx * hx) / (2.0 * sigma * sigma));

/* normalization */
sum = conv[0];
for (i=1; i<=length; i++)
    sum = sum + 2.0 * conv[i];
for (i=0; i<=length; i++)
    conv[i] = conv[i] / sum;

/* allocate storage for a row */
alloc_vector (&help, nx+length+length);

for (j=1; j<=ny; j++)
    {
    /* copy in row vector */
    for (i=1; i<=nx; i++)
        help[i+length-1] = f[i][j];

    /* assign boundary conditions */
    if (bc == 0) /* Dirichlet boundary conditions */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = 0.0;
           help[nx+length-1+p] = 0.0;
           }
    else if (bc == 1) /* reflecting b.c. */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = help[length+p-1];
           help[nx+length-1+p] = help[nx+length-p];
           }
    else if (bc == 2) /* periodic b.c. */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = help[nx+length-p];
           help[nx+length-1+p] = help[length+p-1];
           }

    /* convolution step */
    for (i=length; i<=nx+length-1; i++)
        {
        /* calculate convolution */
        sum = conv[0] * help[i];
        for (p=1; p<=length; p++)
            sum = sum + conv[p] * (help[i+p] + help[i-p]);
        /* write back */
        f[i-length+1][j] = sum;
        }
    } /* for j */

/* disallocate storage for a row */
disalloc_vector (help, nx+length+length);

/* disallocate convolution vector */
disalloc_vector (conv, length + 1);

/* ------------------------ diffusion in y direction -------------------- */

/* calculate length of convolution vector */
length = (long)(precision * sigma / hy) + 1;
if ((bc != 0) && (length > ny))
   {
   printf("gauss_conv: sigma too large \n"); 
   exit(0);
   }

/* allocate storage for convolution vector */
alloc_vector (&conv, length + 1);

/* calculate entries of convolution vector */
for (j=0; j<=length; j++)
    conv[j] = 1 / (sigma * sqrt(2.0 * 3.1415927)) 
              * exp (- (j * j * hy * hy) / (2.0 * sigma * sigma));

/* normalization */
sum = conv[0];
for (j=1; j<=length; j++)
    sum = sum + 2.0 * conv[j];
for (j=0; j<=length; j++)
    conv[j] = conv[j] / sum;

/* allocate storage for a row */
alloc_vector (&help, ny+length+length);

for (i=1; i<=nx; i++)
    {
    /* copy in column vector */
    for (j=1; j<=ny; j++)
        help[j+length-1] = f[i][j];

    /* assign boundary conditions */
    if (bc == 0) /* Dirichlet boundary conditions */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = 0.0;
           help[ny+length-1+p] = 0.0;
           }
    else if (bc == 1) /* reflecting b.c. */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = help[length+p-1];
           help[ny+length-1+p] = help[ny+length-p];
           }
    else if (bc == 2) /* periodic b.c. */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = help[ny+length-p];
           help[ny+length-1+p] = help[length+p-1];
           } 
 
    /* convolution step */
    for (j=length; j<=ny+length-1; j++)
        {
        /* calculate convolution */
        sum = conv[0] * help[j];
        for (p=1; p<=length; p++)
            sum = sum + conv[p] * (help[j+p] + help[j-p]);
        /* write back */
        f[i][j-length+1] = sum;
        }
    } /* for i */

/* disallocate storage for a row */
disalloc_vector (help, ny+length+length);

/* disallocate convolution vector */
disalloc_vector (conv, length+1);

return;

} /* gauss_conv */

/* ------------------------------------------------------------------------ */

void struct_tensor

     (float    ***v,      /* image !! gets smoothed on exit !! */
      long     nc,        /* number of channels */
      long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      float    hx,        /* pixel size in x direction */
      float    hy,        /* pixel size in y direction */
      float    sigma,     /* noise scale */
      float    rho,       /* integration scale */
      float    **dxx,     /* element of structure tensor, output */
      float    **dxy,     /* element of structure tensor, output */
      float    **dyy)     /* element of structure tensor, output */

/*
 Calculates the structure tensor.
*/

{
long    i, j, m;              /* loop variables */
float   dv_dx, dv_dy;         /* derivatives of v */
float   alpha;                /* smoothing time for 3*3 derivative mask */
float   cw;                   /* channel weight */
float   h1, h2, h3, h4;       /* time saver */


/* ---- initialization ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     dxx[i][j] = 0.0;
     dxy[i][j] = 0.0;
     dyy[i][j] = 0.0;
     }

alpha = 0.2;         /* has to be between 0 and 1/3 */
cw    = 1.0 / nc;    /* equal weights to all channels */
h1    = cw * (1.0 - 2.0 * alpha) / (2.0 * hx);
h2    = cw * alpha / (2.0 * hx);
h3    = cw * (1.0 - 2.0 * alpha) / (2.0 * hy);
h4    = cw * alpha / (2.0 * hy);

for (m=0; m<=nc-1; m++)
    {
    /* ---- smoothing at noise scale ---- */

    if (sigma > 0.0)
       gauss_conv (sigma, nx, ny, hx, hy, 3.0, 1, v[m]);  /* refl. b.c. */


    /* ---- building tensor product ---- */

    dummies (v[m], nx, ny);

    for (i=1; i<=nx; i++)
     for (j=1; j<=ny; j++)
         {
         /* modified Sobel derivatives */
         dv_dx =   h1 * ( v[m][i+1][j]   - v[m][i-1][j] )
                 + h2 * ( v[m][i+1][j-1] - v[m][i-1][j-1]
                        + v[m][i+1][j+1] - v[m][i-1][j+1] );
         dv_dy =   h3 * ( v[m][i][j+1]   - v[m][i][j-1] )
                 + h4 * ( v[m][i-1][j+1] - v[m][i-1][j-1]
                        + v[m][i+1][j+1] - v[m][i+1][j-1] );
         /* tensor product */
         dxx[i][j] = dxx[i][j] + dv_dx * dv_dx;
         dxy[i][j] = dxy[i][j] + dv_dx * dv_dy;
         dyy[i][j] = dyy[i][j] + dv_dy * dv_dy;
         }
    } /* for m */


/* ---- smoothing at integration scale, Dirichlet b.c. ---- */

if (rho > 0.0)
   {
   gauss_conv (rho, nx, ny, hx, hy, 3.0, 0, dxx);
   gauss_conv (rho, nx, ny, hx, hy, 3.0, 0, dxy);
   gauss_conv (rho, nx, ny, hx, hy, 3.0, 0, dyy);
   }

return;

} /* struct_tensor */

/* ------------------------------------------------------------------------ */

void PA_trans

     (float a11,        /* coeffs of (2*2)-matrix */
      float a12,        /* coeffs of (2*2)-matrix */
      float a22,        /* coeffs of (2*2)-matrix */
      float *c,         /* 1. comp. of 1. eigenvector, output */
      float *s,         /* 2. comp. of 1. eigenvector, output */
      float *lam1,      /* larger  eigenvalue, output */
      float *lam2)      /* smaller eigenvalue, output */

/*
 Principal axis transformation, checked for correctness.
*/

{
float  help, norm;    /* time savers */


/* ---- compute eigenvalues and eigenvectors ---- */

help  = sqrt (powf (a22-a11, 2.0) + 4.0 * a12 * a12);

if (help == 0.0)
   /* isotropic situation, eigenvectors arbitrary */
   {
   *lam1 = *lam2 = a11;
   *c = 1.0;
   *s = 0.0;
   }

else if (a11 > a22)
   {
   *lam1 = 0.5 * (a11 + a22 + help);
   *lam2 = 0.5 * (a11 + a22 - help);
   *c = a11 - a22 + help;
   *s = 2.0 * a12;
   }

else
   {
   *lam1 = 0.5 * (a11 + a22 + help);
   *lam2 = 0.5 * (a11 + a22 - help);
   *c = 2.0 * a12;
   *s = a22 - a11 + help;
   }


/* ---- normalize eigenvectors ---- */

norm = sqrt (*c * *c + *s * *s);
if (norm > 0.0)
   {
   *c = *c / norm;
   *s = *s / norm;
   }
else
   {
   *c = 1.0;
   *s = 0.0;
   }

return;

} /* PA_trans */

/*--------------------------------------------------------------------------*/

void corner_detection 

     (float    ***f,      /* image, unaltered */
      long     nc,        /* number of channels of u*/
      long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      float    hx,        /* pixel size in x direction */
      float    hy,        /* pixel size in y direction */
      long     type,      /* type of corner detector */
      float    T,         /* threshold */
      float    sigma,     /* noise scale */
      float    rho,       /* integration scale */
      float    **v)       /* corner locations */

/*
 calculates structure tensor based corner detector; 
 v[i][j]=255 for corners; 
 v[i][j]=0   else; 
*/

{
long    i, j, m;               /* loop variables */
float   ***u;                  /* loop */
float   **dxx, **dxy, **dyy;   /* tensor components */
float   **w;                   /* field for corner detection */
float   c, s;                  /* cosine, sine */
float   lam1, lam2;            /* eigenvalues */
float   det, trace;            /* determinant, trace */


/* allocate storage */
alloc_cubix  (&u, nc, nx+2, ny+2);
alloc_matrix (&dxx, nx+2, ny+2);
alloc_matrix (&dxy, nx+2, ny+2);
alloc_matrix (&dyy, nx+2, ny+2);
alloc_matrix (&w, nx+2, ny+2);

/* copy f into u */
for (m=0; m<=nc-1; m++)
 for (i=1; i<=nx; i++)
  for (j=1; j<=ny; j++)
      u[m][i][j] = f[m][i][j] ;

/* calculate structure tensor of u */
struct_tensor (u, nc, nx, ny, hx, hy, sigma, rho, dxx, dxy, dyy);

if (type == 0)
   /* Rohr: det(J) */
   for (i=1; i<=nx; i++)
    for (j=1; j<=ny; j++)
        {
        w[i][j] = dxx[i][j] * dyy[i][j] - dxy[i][j] * dxy[i][j];        
        if (w[i][j] <= T)
           w[i][j] = 0.0;
        } 

if (type == 1)
   /* Tomasi-Kanade: smaller eigenvalue */
   for (i=1; i<=nx; i++)
    for (j=1; j<=ny; j++)
        {
        PA_trans (dxx[i][j], dxy[i][j], dyy[i][j], &c, &s, &lam1, &lam2);
        w[i][j] = lam2;
        if (w[i][j] <= T)
           w[i][j] = 0.0;
        }

if (type == 2)
   /* Foerstner-Harris: det(J)/trace(J) */
   for (i=1; i<=nx; i++)
    for (j=1; j<=ny; j++)
        {
        trace = dxx[i][j] + dyy[i][j];
        det   = dxx[i][j] * dyy[i][j] - dxy[i][j] * dxy[i][j];
        if (trace > T) 
           w[i][j] = det / trace;
        else
           w[i][j] = 0.0;
        }

/* search for maximum */
dummies (w, nx, ny);
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     if ((w[i][j] > w[i+1][j]) && (w[i][j] > w[i-1][j]) && 
         (w[i][j] > w[i][j+1]) && (w[i][j] > w[i][j-1]) &&
         (w[i+1][j] * w[i-1][j] * w[i][j+1] * w[i][j-1] != 0.0) )
        v[i][j] = 255.0;
     else
        v[i][j] = 0.0;


/* free storage */
disalloc_cubix (u, nc, nx+2, ny+2);
disalloc_matrix (dxx, nx+2, ny+2);
disalloc_matrix (dxy, nx+2, ny+2);
disalloc_matrix (dyy, nx+2, ny+2);
disalloc_matrix (w,   nx+2, ny+2);

return;

} /* corner_detection */

/*--------------------------------------------------------------------------*/

void draw_corners 

     (float    ***u,      /* image, altered ! */
      long     nc,        /* number of channels of u*/
      long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      float    **v)       /* image with corner locations */

/*
 draws corners into the image u, at locations specified by v 
 v[i][j]=255 for corners; 
 v[i][j]=0   else; 
*/

{
long    i, j, m;              /* loop variables */
long    k, l;                 /* loop variables */

for (m=0; m<=nc-1; m++)
    dummies (u[m], nx, ny);

for (m=0; m<=nc-1; m++)
 for (i=5; i<=nx-4; i++)
  for (j=5; j<=ny-4; j++)
      if (v[i][j] == 255.0)
         /* draw corner */
         for (k=i-4; k<=i+4; k++)
           for (l=j-4; l<=j+4; l++)
               { 
               /* black outer circle */
               if ((k-i)*(k-i)+(l-j)*(l-j) <= 20) 
                  u[m][k][l] = 0.0; 
               /* white interior */
               if ((k-i)*(k-i)+(l-j)*(l-j) <= 6) 
                  u[m][k][l] = 255.0; 
               }

return;

} /* draw_corners */

/*--------------------------------------------------------------------------*/

int main ()

{
char   row[80];              /* for reading data */
char   in[80];               /* for reading data */
char   out[80];              /* for reading data */
float  ***u;                 /* input image */
float  **v;                  /* image with corner location */
long   i, j, m;              /* loop variables */
long   nx, ny;               /* image size in x, y direction */
long   nc;                   /* number of channels */
long   type;                 /* type of corner detector */ 
long   count;                /* number of corners */
float  T;                    /* threshold */
float  sigma;                /* noise scale */
float  rho;                  /* integration scale */
float  max, min;             /* largest, smallest grey value */
unsigned char byte;          /* for data conversion */
FILE   *inimage, *outimage;  /* input file, output file */


printf("\n");
printf("CORNER DETECTION WITH THE STRUCTURE TENSOR\n\n");

/* ---- read input image (pgm format P5 or ppm format P6) ---- */

/* read image name */
printf("input image (pgm, ppm):           ");
gets (in); 

/* open pgm file and read header */
inimage = fopen(in,"r");
fgets (row, 300, inimage);                       /* image type: P5 or P6 */
if ((row[0]=='P') && (row[1]=='5'))
   nc = 1;                                       /* P5: grey scale image */
else if ((row[0]=='P') && (row[1]=='6'))
   nc = 3;                                       /* P6: colour image */
else
   {
   printf ("unknown image format");
   exit(0);
   }
fgets (row, 300, inimage);                       /* first comment row */
while (row[0]=='#') fgets(row, 300, inimage);    /* further comment rows */
sscanf (row, "%ld %ld", &nx, &ny);               /* image size */
fgets (row, 300, inimage);                       /* # grey values */

/* allocate storage */
alloc_cubix (&u, nc, nx+2, ny+2);

/* read image data */
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
  for (m=0; m<=nc-1; m++)
      u[m][i][j] = (float) getc (inimage);
fclose(inimage);


/* ---- read other parameters ---- */

printf("corner detector:\n");
printf("  (0) Rohr:             det(J)\n");
printf("  (1) Tomasi-Kanade:    lambda_2\n");
printf("  (2) Foerstner-Harris: det(J)/tr(J)\n");
printf("your choice:                      ");
gets(row);  sscanf(row, "%ld", &type);
printf("threshold T (>=0):                ");
gets(row);  sscanf(row, "%f", &T);
printf("noise scale sigma (>=0):          ");
gets(row);  sscanf(row, "%f", &sigma);
printf("integration scale rho (>=0):      ");
gets(row);  sscanf(row, "%f", &rho);
if (nc == 1) 
   printf("output image (pgm):               ");
else if (nc == 3) 
   printf("output image (ppm):               ");
gets(out);
printf("\n");


/* ---- process image ---- */

/* allocate storage for corner detector */
alloc_matrix (&v, nx+2, ny+2);

/* corner detection */
corner_detection (u, nc, nx, ny, 1.0, 1.0, type, T, sigma, rho, v); 

/* insert corner marks into the image */
draw_corners (u, nc, nx, ny, v);

/* count corners */
count = 0;
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     if (v[i][j] == 255.0)
        count = count + 1;
printf("number of corners:   %3ld\n", count);


/* ---- write output image ---- */

/* open file and write header (incl. filter parameters) */
outimage = fopen (out, "w");
if (nc == 1)
   fprintf (outimage, "P5 \n");
else if (nc == 3)
   fprintf (outimage, "P6 \n");
fprintf (outimage, "# structure tensor analysis\n");
if (type == 0) fprintf (outimage, "# Rohr corner detector \n");
if (type == 1) fprintf (outimage, "# Tomasi-Kanade corner detector \n");
if (type == 2) fprintf (outimage, "# Foerstner-Harris corner detector \n");
fprintf (outimage, "# initial image:  %s\n", in);
fprintf (outimage, "# sigma:          %8.4f\n", sigma);
fprintf (outimage, "# rho:            %8.4f\n", rho);
fprintf (outimage, "# T:              %8.4f\n", T);
fprintf (outimage, "# corners:        %8ld\n", count);
fprintf (outimage, "%ld %ld \n255\n", nx, ny);

/* write image data and close file */
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
  for (m=0; m<=nc-1; m++)
     {
     if (u[m][i][j] < 0.0)
        byte = (unsigned char)(0.0);
     else if (u[m][i][j] > 255.0)
        byte = (unsigned char)(255.0);
     else
        byte = (unsigned char)(u[m][i][j]);
     fwrite (&byte, sizeof(unsigned char), 1, outimage);
     }
fclose(outimage);
printf("output image %s successfully written\n\n", out);


/* ---- disallocate storage ---- */

disalloc_cubix (u, nc, nx+2, ny+2);
disalloc_matrix (v, nx+2, ny+2);

return(0);
}
