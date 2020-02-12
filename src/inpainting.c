#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define float double

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                           EED-BASED INPAINTING                           */
/*                                                                          */
/*                   (Copyright Joachim Weickert, 11/2008)                  */
/*                                                                          */
/*--------------------------------------------------------------------------*/

/*
 features:
 - parabolic PDEs
 - new space discretisation for EED
 - semi-implicit scheme with CG solver
*/

/*--------------------------------------------------------------------------*/

void alloc_vector

    (float** vector, /* vector */
        long n)      /* size */

/* allocates storage for a vector of size n */

{
    *vector = (float*)malloc(n * sizeof(float));
    if (*vector == NULL) {
        printf("alloc_vector: not enough storage available\n");
        exit(1);
    }
    return;
}

/*--------------------------------------------------------------------------*/

void alloc_matrix

    (float*** matrix, /* matrix */
        long nx,      /* size in x direction */
        long ny)      /* size in y direction */

/* allocates storage for matrix of size nx * ny */

{
    long i;

    *matrix = (float**)malloc(nx * sizeof(float*));
    if (*matrix == NULL) {
        printf("alloc_matrix: not enough storage available\n");
        exit(1);
    }
    for (i = 0; i < nx; i++) {
        (*matrix)[i] = (float*)malloc(ny * sizeof(float));
        if ((*matrix)[i] == NULL) {
            printf("alloc_matrix: not enough storage available\n");
            exit(1);
        }
    }
    return;
}

/*--------------------------------------------------------------------------*/

void disalloc_vector

    (float* vector, /* vector */
        long n)     /* size */

/* disallocates storage for a vector of size n */

{
    free(vector);
    return;
}

/*--------------------------------------------------------------------------*/

void disalloc_matrix

    (float** matrix, /* matrix */
        long nx,     /* size in x direction */
        long ny)     /* size in y direction */

/* disallocates storage for matrix of size nx * ny */

{
    long i;
    for (i = 0; i < nx; i++)
        free(matrix[i]);
    free(matrix);
    return;
}

/*--------------------------------------------------------------------------*/

void dummies_Dirichlet

    (float** v,  /* image matrix */
        long nx, /* size in x direction */
        long ny) /* size in y direction */

/* creates homogeneous Dirichlet boundaries */

{
    long i, j; /* loop variables */

    for (i = 1; i <= nx; i++) {
        v[i][0] = 0.0;
        v[i][ny + 1] = 0.0;
    }

    for (j = 0; j <= ny + 1; j++) {
        v[0][j] = 0.0;
        v[nx + 1][j] = 0.0;
    }
    return;
}

/*--------------------------------------------------------------------------*/

void dummies_Neumann

    (float** v,  /* image matrix */
        long nx, /* size in x direction */
        long ny) /* size in y direction */

/* creates homogeneous Neumann boundaries by mirroring */

{
    long i, j; /* loop variables */

    for (i = 1; i <= nx; i++) {
        v[i][0] = v[i][1];
        v[i][ny + 1] = v[i][ny];
    }

    for (j = 0; j <= ny + 1; j++) {
        v[0][j] = v[1][j];
        v[nx + 1][j] = v[nx][j];
    }
    return;
}

/* ---------------------------------------------------------------------- */

void gauss_conv

    (float sigma,        /* standard deviation of Gaussian */
        long nx,         /* image dimension in x direction */
        long ny,         /* image dimension in y direction */
        float hx,        /* pixel size in x direction */
        float hy,        /* pixel size in y direction */
        float precision, /* cutoff at precision * sigma */
        long bc,         /* type of boundary condition */
                         /* 0=Dirichlet, 1=reflecing, 2=periodic */
        float** f)       /* input: original image ;  output: smoothed */

/*
 Gaussian convolution.
*/

{
    long i, j, p; /* loop variables */
    long length;  /* convolution vector: 0..length */
    float sum;    /* for summing up */
    float* conv;  /* convolution vector */
    float* help;  /* row or column with dummy boundaries */

    /* ------------------------ convolution in x direction -------------------- */

    /* calculate length of convolution vector */
    length = (long)(precision * sigma / hx) + 1;
    if ((bc != 0) && (length > nx)) {
        printf("gauss_conv: sigma too large \n");
        exit(0);
    }

    /* allocate storage for convolution vector */
    alloc_vector(&conv, length + 1);

    /* calculate entries of convolution vector */
    for (i = 0; i <= length; i++)
        conv[i] = 1 / (sigma * sqrt(2.0 * 3.1415926))
            * exp(-(i * i * hx * hx) / (2.0 * sigma * sigma));

    /* normalisation */
    sum = conv[0];
    for (i = 1; i <= length; i++)
        sum = sum + 2.0 * conv[i];
    for (i = 0; i <= length; i++)
        conv[i] = conv[i] / sum;

    /* allocate storage for a row */
    alloc_vector(&help, nx + length + length);

    for (j = 1; j <= ny; j++) {
        /* copy in row vector */
        for (i = 1; i <= nx; i++)
            help[i + length - 1] = f[i][j];

        /* assign boundary conditions */
        if (bc == 0) /* Dirichlet boundary conditions */
            for (p = 1; p <= length; p++) {
                help[length - p] = 0.0;
                help[nx + length - 1 + p] = 0.0;
            }
        else if (bc == 1) /* reflecting b.c. */
            for (p = 1; p <= length; p++) {
                help[length - p] = help[length + p - 1];
                help[nx + length - 1 + p] = help[nx + length - p];
            }
        else if (bc == 2) /* periodic b.c. */
            for (p = 1; p <= length; p++) {
                help[length - p] = help[nx + length - p];
                help[nx + length - 1 + p] = help[length + p - 1];
            }

        /* convolution step */
        for (i = length; i <= nx + length - 1; i++) {
            /* calculate convolution */
            sum = conv[0] * help[i];
            for (p = 1; p <= length; p++)
                sum = sum + conv[p] * (help[i + p] + help[i - p]);
            /* write back */
            f[i - length + 1][j] = sum;
        }
    } /* for j */

    /* disallocate storage for a row */
    disalloc_vector(help, nx + length + length);

    /* disallocate convolution vector */
    disalloc_vector(conv, length + 1);

    /* ------------------------ convolution in y direction -------------------- */

    /* calculate length of convolution vector */
    length = (long)(precision * sigma / hy) + 1;
    if ((bc != 0) && (length > ny)) {
        printf("gauss_conv: sigma too large \n");
        exit(0);
    }

    /* allocate storage for convolution vector */
    alloc_vector(&conv, length + 1);

    /* calculate entries of convolution vector */
    for (j = 0; j <= length; j++)
        conv[j] = 1 / (sigma * sqrt(2.0 * 3.1415927))
            * exp(-(j * j * hy * hy) / (2.0 * sigma * sigma));

    /* normalisation */
    sum = conv[0];
    for (j = 1; j <= length; j++)
        sum = sum + 2.0 * conv[j];
    for (j = 0; j <= length; j++)
        conv[j] = conv[j] / sum;

    /* allocate storage for a row */
    alloc_vector(&help, ny + length + length);

    for (i = 1; i <= nx; i++) {
        /* copy in column vector */
        for (j = 1; j <= ny; j++)
            help[j + length - 1] = f[i][j];

        /* assign boundary conditions */
        if (bc == 0) /* Dirichlet boundary conditions */
            for (p = 1; p <= length; p++) {
                help[length - p] = 0.0;
                help[ny + length - 1 + p] = 0.0;
            }
        else if (bc == 1) /* reflecting b.c. */
            for (p = 1; p <= length; p++) {
                help[length - p] = help[length + p - 1];
                help[ny + length - 1 + p] = help[ny + length - p];
            }
        else if (bc == 2) /* periodic b.c. */
            for (p = 1; p <= length; p++) {
                help[length - p] = help[ny + length - p];
                help[ny + length - 1 + p] = help[length + p - 1];
            }

        /* convolution step */
        for (j = length; j <= ny + length - 1; j++) {
            /* calculate convolution */
            sum = conv[0] * help[j];
            for (p = 1; p <= length; p++)
                sum = sum + conv[p] * (help[j + p] + help[j - p]);
            /* write back */
            f[i][j - length + 1] = sum;
        }
    } /* for i */

    /* disallocate storage for a row */
    disalloc_vector(help, ny + length + length);

    /* disallocate convolution vector */
    disalloc_vector(conv, length + 1);

    return;

} /* gauss_conv */

/* ------------------------------------------------------------------------ */

void struct_tensor

    (float** v,      /* image !! gets smoothed on exit !! */
        long nx,     /* image dimension in x direction */
        long ny,     /* image dimension in y direction */
        float hx,    /* pixel size in x direction */
        float hy,    /* pixel size in y direction */
        float sigma, /* noise scale */
        float rho,   /* integration scale */
        float** dxx, /* element of structure tensor, output */
        float** dxy, /* element of structure tensor, output */
        float** dyy) /* element of structure tensor, output */

/*
 Computes the structure tensor at intermediate grid points
 [i+1/2][j+1/2] for i=1,...,nx-1 and j=1,...,ny-1.
*/

{
    long i, j;                /* loop variables */
    float vxo, vxp, vyo, vyp; /* derivatives of v */
    float vxc, vyc;           /* averaged derivatives of v */
    float w1, w2;             /* time savers */

    /* ---- smoothing at noise scale, reflecting b.c. ---- */

    if (sigma > 0.0)
        gauss_conv(sigma, nx, ny, hx, hy, 5.0, 1, v);

    /* ---- building tensor product ---- */

    w1 = 1.0 / hx;
    w2 = 1.0 / hy;
    dummies_Neumann(v, nx, ny);

    for (i = 1; i <= nx - 1; i++)
        for (j = 1; j <= ny - 1; j++) {
            vxo = w1 * (v[i + 1][j] - v[i][j]);
            vxp = w1 * (v[i + 1][j + 1] - v[i][j + 1]);
            vyo = w2 * (v[i][j + 1] - v[i][j]);
            vyp = w2 * (v[i + 1][j + 1] - v[i + 1][j]);
            vxc = 0.5 * (vxo + vxp);
            vyc = 0.5 * (vyo + vyp);
            dxx[i][j] = 0.25 * (vxo * vxo + vxp * vxp) + 0.5 * (vxc * vxc);
            dyy[i][j] = 0.25 * (vyo * vyo + vyp * vyp) + 0.5 * (vyc * vyc);
            dxy[i][j] = vxc * vyc;
        }

    /* ---- smoothing at integration scale, Dirichlet b.c. ---- */

    if (rho > 0.0) {
        gauss_conv(rho, nx - 1, ny - 1, hx, hy, 5.0, 0, dxx);
        gauss_conv(rho, nx - 1, ny - 1, hx, hy, 5.0, 0, dxy);
        gauss_conv(rho, nx - 1, ny - 1, hx, hy, 5.0, 0, dyy);
    }

    return;

} /* struct_tensor */

/* ------------------------------------------------------------------------ */

void PA_trans

    (float a11,      /* coeffs of (2*2)-matrix */
        float a12,   /* coeffs of (2*2)-matrix */
        float a22,   /* coeffs of (2*2)-matrix */
        float* c,    /* 1. comp. of 1. eigenvector, output */
        float* s,    /* 2. comp. of 1. eigenvector, output */
        float* lam1, /* larger  eigenvalue, output */
        float* lam2) /* smaller eigenvalue, output */

/*
 Principal axis transformation, checked for correctness.
*/

{
    float help, norm; /* time savers */

    /* ---- compute eigenvalues and eigenvectors ---- */

    help = sqrt(pow(a22 - a11, 2.0) + 4.0 * a12 * a12);

    if (help == 0.0)
    /* isotropic situation, eigenvectors arbitrary */
    {
        *lam1 = *lam2 = a11;
        *c = 1.0;
        *s = 0.0;
    }

    else if (a11 > a22) {
        *lam1 = 0.5 * (a11 + a22 + help);
        *lam2 = 0.5 * (a11 + a22 - help);
        *c = a11 - a22 + help;
        *s = 2.0 * a12;
    }

    else {
        *lam1 = 0.5 * (a11 + a22 + help);
        *lam2 = 0.5 * (a11 + a22 - help);
        *c = 2.0 * a12;
        *s = a22 - a11 + help;
    }

    /* ---- normalize eigenvectors ---- */

    norm = sqrt(*c * *c + *s * *s);
    if (norm >= 0.000001) {
        *c = *c / norm;
        *s = *s / norm;
    } else {
        *c = 1.0;
        *s = 0.0;
    }

    return;

} /* PA_trans */

/* ------------------------------------------------------------------------ */

void PA_backtrans

    (float c,       /* 1. comp. of 1. eigenvector */
        float s,    /* 2. comp. of 1. eigenvector */
        float lam1, /* 1. eigenvalue */
        float lam2, /* 2. eigenvalue */
        float* a11, /* coeff. of (2*2)-matrix, output */
        float* a12, /* coeff. of (2*2)-matrix, output */
        float* a22) /* coeff. of (2*2)-matrix, output */

/*
 Principal axis backtransformation of a symmetric (2*2)-matrix.
 A = U * diag(lam1, lam2) * U_transpose with U = (v1 | v2)
 v1 = (c, s) is first eigenvector
*/

{

    *a11 = c * c * lam1 + s * s * lam2;
    *a22 = lam1 + lam2 - *a11; /* trace invariance */
    *a12 = c * s * (lam1 - lam2);

    return;

} /* PA_backtrans */

/*--------------------------------------------------------------------------*/

void diff_tensor

    (long dtype,      /* type of diffusivity */
        float lambda, /* contrast parameter */
        long nx,      /* image dimension in x direction */
        long ny,      /* image dimension in y direction */
        float** dxx,  /* in: structure tensor el., out: diff. tensor el. */
        float** dxy,  /* in: structure tensor el., out: diff. tensor el. */
        float** dyy)  /* in: structure tensor el., out: diff. tensor el. */

/*
 Calculates the diffusion tensor by means of the structure tensor.
*/

{
    long i, j;        /* loop variables */
    float help;       /* time saver */
    float c, s;       /* specify first eigenvector */
    float mu1, mu2;   /* eigenvalues of structure tensor */
    float lam1, lam2; /* eigenvalues of diffusion tensor */

    /* ---- compute diffusion tensor ---- */

    help = 1.0 / (lambda * lambda);

    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++) {
            /* principal axis transformation */
            PA_trans(dxx[i][j], dxy[i][j], dyy[i][j], &c, &s, &mu1, &mu2);

            /* calculate eigenvalues */
            if (dtype == 0)
                /* Charbonnier diffusivity */
                lam1 = 1.0 / sqrt(1.0 + mu1 * help);
            else if (dtype == 1)
                /* Weickert diffusivity */
                lam1 = 1.0 - exp(-3.31488 / powf(mu1 * help, 4.0));
            lam2 = 1.0;

            /* principal axis backtransformation */
            PA_backtrans(c, s, lam1, lam2, &dxx[i][j], &dxy[i][j], &dyy[i][j]);
        }

    /* ---- assign dummy boundaries (no flux) ---- */

    dummies_Dirichlet(dxx, nx, ny);
    dummies_Dirichlet(dxy, nx, ny);
    dummies_Dirichlet(dyy, nx, ny);

    return;

} /* diff_tensor */

/* -------------------------------------------------------------------------*/

float sgn

    (float x) /* argument */

/*
 sign function
*/

{
    float sign; /* auxiliary variable */

    if (x > 0)
        sign = 1.0;
    else if (x < 0)
        sign = -1.0;
    else
        sign = 0.0;

    return (sign);
}

/*--------------------------------------------------------------------------*/

void weights

    (float** dxx,    /* entry [1,1] of structure tensor, unchanged */
        float** dxy, /* entry [1,2] of structure tensor, unchanged */
        float** dyy, /* entry [2,2] of structure tensor, unchanged */
        float** a,   /* confidence map, unchanged */
        long nx,     /* image dimension in x direction */
        long ny,     /* image dimension in y direction */
        float hx,    /* pixel size in x direction */
        float hy,    /* pixel size in y direction */
        float alpha, /* dissipativity parameter, in [0,0.5] */
        float gamma, /* nonnegativity parameter, in [0,1] */
        float** woo, /* weights in [i,j], output */
        float** wpo, /* weights in [i+1,j], output */
        float** wmo, /* weights in [i-1,j], output */
        float** wop, /* weights in [i,j+1], output */
        float** wom, /* weights in [i,j-1], output */
        float** wpp, /* weights in [i+1,j+1], output */
        float** wmm, /* weights in [i-1,j-1], output */
        float** wpm, /* weights in [i+1,j-1], output */
        float** wmp) /* weights in [i-1,j+1], output */

/*
 computes weights arising from the discrete divergence expression
*/

{
    long i, j;           /* loop variables */
    float** beta;        /* space-variant numerical parameter */
    float rxx, ryy, rxy; /* time savers */

    /* ---- allocate storage ---- */

    alloc_matrix(&beta, nx + 1, ny + 1);

    /* ---- initialisations ---- */

    rxx = 1.0 / (2.0 * hx * hx);
    ryy = 1.0 / (2.0 * hy * hy);
    rxy = 1.0 / (2.0 * hx * hy);

    for (i = 0; i <= nx; i++)
        for (j = 0; j <= ny; j++)
            beta[i][j] = gamma * (1.0 - 2.0 * alpha) * sgn(dxy[i][j]);

    for (i = 0; i <= nx + 1; i++)
        for (j = 0; j <= ny + 1; j++) {
            woo[i][j] = 0.0;
            wpo[i][j] = wop[i][j] = wpp[i][j] = wpm[i][j] = 0.0;
            wmo[i][j] = wom[i][j] = wmm[i][j] = wmp[i][j] = 0.0;
        }

    /* ---- weights ---- */

    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++) {
            wpo[i][j] = (1.0 - alpha) * rxx * (dxx[i][j] + dxx[i][j - 1])
                - alpha * ryy * (dyy[i][j] + dyy[i][j - 1])
                - rxy * (beta[i][j] * dxy[i][j] + beta[i][j - 1] * dxy[i][j - 1]);

            wop[i][j] = (1.0 - alpha) * ryy * (dyy[i][j] + dyy[i - 1][j])
                - alpha * rxx * (dxx[i][j] + dxx[i - 1][j])
                - rxy * (beta[i][j] * dxy[i][j] + beta[i - 1][j] * dxy[i - 1][j]);

            wpp[i][j] = alpha * rxx * dxx[i][j] + alpha * ryy * dyy[i][j]
                + (beta[i][j] + 1.0) * rxy * dxy[i][j];

            wpm[i][j] = alpha * rxx * dxx[i][j - 1] + alpha * ryy * dyy[i][j - 1]
                + (beta[i][j - 1] - 1.0) * rxy * dxy[i][j - 1];
        }

    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++) {
            wmo[i][j] = wpo[i - 1][j];
            wom[i][j] = wop[i][j - 1];
            wmm[i][j] = wpp[i - 1][j - 1];
            wmp[i][j] = wpm[i - 1][j + 1];

            woo[i][j] = -wpo[i][j] - wop[i][j] - wpp[i][j] - wpm[i][j] - wmo[i][j] - wom[i][j]
                - wmm[i][j] - wmp[i][j];
        }

    /* don't change anything at interpolation points */
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
            if (a[i][j] > 0.00001) {
                woo[i][j] = 0.0;
                wpo[i][j] = wop[i][j] = wpp[i][j] = wpm[i][j] = 0.0;
                wmo[i][j] = wom[i][j] = wmm[i][j] = wmp[i][j] = 0.0;
            }

    /* ---- disallocate storage ---- */

    disalloc_matrix(beta, nx + 1, ny + 1);

    return;

} /* weights */

/*--------------------------------------------------------------------------*/

void weights_to_matrix

    (long nx,        /* image dimension in x direction */
        long ny,     /* image dimension in y direction */
        float ht,    /* time step size */
        float** boo, /* entries for [i,j],     changed */
        float** bpo, /* entries for [i+1,j],   changed */
        float** bmo, /* entries for [i-1,j],   changed */
        float** bop, /* entries for [i,j+1],   changed */
        float** bom, /* entries for [i,j-1],   changed */
        float** bpp, /* entries for [i+1,j+1], changed */
        float** bmm, /* entries for [i-1,j-1], changed */
        float** bpm, /* entries for [i+1,j-1], changed */
        float** bmp) /* entries for [i-1,j+1], changed */

/*
 input:  symmetric, nonadiagonal matrix A with diagonal boo and
         off-diagonals bpo, bmo, bop, bom, bpp, bmm, bpm, bmp
         corresponding to the 8 neighbours
 output: symmetric, nonadiagonal matrix B = I - ht * A with diagonal
         boo and 8 off-diagonals bpo, bmo, bop, bom, bpp, bmm, bpm, bmp.
*/

{
    long i, j; /* loop variables */

    for (i = 0; i <= nx + 1; i++)
        for (j = 0; j <= ny + 1; j++) {
            bpo[i][j] = -ht * bpo[i][j];
            bmo[i][j] = -ht * bmo[i][j];
            bop[i][j] = -ht * bop[i][j];
            bom[i][j] = -ht * bom[i][j];
            bpp[i][j] = -ht * bpp[i][j];
            bmm[i][j] = -ht * bmm[i][j];
            bpm[i][j] = -ht * bpm[i][j];
            bmp[i][j] = -ht * bmp[i][j];
            boo[i][j] = 1.0 - ht * boo[i][j];
        }

    return;

} /* weights_to_matrix */

/*--------------------------------------------------------------------------*/

float norm

    (long nx,      /* image dimension in x direction */
        long ny,   /* image dimension in y direction */
        float** u) /* image, unchanged */

/*
  computes the normalised L2 norm of the image u
*/

{
    long i, j; /* loop variables */
    float aux; /* auxiliary variable */

    /* ---- sum up squared contributions ---- */

    aux = 0.0;
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
            aux = aux + u[i][j] * u[i][j];

    /* ---- take square root and normalise ---- */

    aux = sqrt(aux) / (nx * ny);

    return (aux);
}

/*--------------------------------------------------------------------------*/

float residue

    (long nx,        /* image dimension in x direction */
        long ny,     /* image dimension in y direction */
        float** boo, /* diagonal entries for [i,j], unchanged */
        float** bpo, /* neighbour entries for [i+1,j], unchanged */
        float** bmo, /* neighbour entries for [i-1,j], unchanged */
        float** bop, /* neighbour entries for [i,j+1], unchanged */
        float** bom, /* neighbour entries for [i,j-1], unchanged */
        float** bpp, /* neighbour entries for [i+1,j+1], unchanged */
        float** bmm, /* neighbour entries for [i-1,j-1], unchanged */
        float** bpm, /* neighbour entries for [i+1,j-1], unchanged */
        float** bmp, /* neighbour entries for [i-1,j+1], unchanged */
        float** f,   /* right hand side, unchanged */
        float** u)   /* approximation to the solution, unchanged */

/*
  computes normalized L^2 norm of the residue of a linear system
  B u = f with a symmetric, nonadiagonal system matrix B that involves
  all eight 2D neighbours.
*/

{
    long i, j;  /* loop variables */
    float r;    /* normalized L2 norm of the residue */
    float help; /* auxiliary variable */

    /* ---- initialisations ---- */

    dummies_Neumann(u, nx, ny);
    r = 0.0;

    /* ---- compute squared L^2 norm of the residue ---- */

    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++) {
            help = f[i][j]
                - (boo[i][j] * u[i][j] + bpo[i][j] * u[i + 1][j] + bmo[i][j] * u[i - 1][j]
                    + bop[i][j] * u[i][j + 1] + bom[i][j] * u[i][j - 1]
                    + bpp[i][j] * u[i + 1][j + 1] + bmm[i][j] * u[i - 1][j - 1]
                    + bpm[i][j] * u[i + 1][j - 1] + bmp[i][j] * u[i - 1][j + 1]);
            r = r + help * help;
        }

    /* ---- take square root and normalise ---- */

    r = sqrt(r) / (nx * ny);

    return (r);
}

/*--------------------------------------------------------------------------*/

float inner_product

    (long nx,      /* image dimension in x direction */
        long ny,   /* image dimension in y direction */
        float** u, /* image 1, unchanged */
        float** v) /* image 2, unchanged */

/*
  computes the inner product of two vectors u and v
*/

{
    long i, j; /* loop variables */
    float aux; /* auxiliary variable */

    aux = 0.0;
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
            aux = aux + u[i][j] * v[i][j];

    return (aux);
}

/*--------------------------------------------------------------------------*/

void matrix_times_vector

    (long nx,        /* image dimension in x direction */
        long ny,     /* image dimension in y direction */
        float** boo, /* matrix diagonal entries for [i,j], unchanged */
        float** bpo, /* neighbour entries for [i+1,j], unchanged */
        float** bmo, /* neighbour entries for [i-1,j], unchanged */
        float** bop, /* neighbour entries for [i,j+1], unchanged */
        float** bom, /* neighbour entries for [i,j-1], unchanged */
        float** bpp, /* neighbour entries for [i+1,j+1], unchanged */
        float** bmm, /* neighbour entries for [i-1,j-1], unchanged */
        float** bpm, /* neighbour entries for [i+1,j-1], unchanged */
        float** bmp, /* neighbour entries for [i-1,j+1], unchanged */
        float** f,   /* vector, unchanged */
        float** u)   /* result, changed */

/*
  computes the product of a symmetric nonadiagonal matrix specified by
  the diagonal boo and the off-diagonals bpo,..., bmp and a vector f
*/

{
    long i, j; /* loop variables */

    dummies_Neumann(f, nx, ny);

    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
            u[i][j] = (boo[i][j] * f[i][j] + bpo[i][j] * f[i + 1][j] + bmo[i][j] * f[i - 1][j]
                + bop[i][j] * f[i][j + 1] + bom[i][j] * f[i][j - 1] + bpp[i][j] * f[i + 1][j + 1]
                + bmm[i][j] * f[i - 1][j - 1] + bpm[i][j] * f[i + 1][j - 1]
                + bmp[i][j] * f[i - 1][j + 1]);

    return;
}

/*--------------------------------------------------------------------------*/

void CG

    (long kmax,      /* max. number of iterations */
        long nx,     /* image dimension in x direction */
        long ny,     /* image dimension in y direction */
        float** boo, /* diagonal entries for [i,j], unchanged */
        float** bpo, /* neighbour entries for [i+1,j], unchanged */
        float** bmo, /* neighbour entries for [i-1,j], unchanged */
        float** bop, /* neighbour entries for [i,j+1], unchanged */
        float** bom, /* neighbour entries for [i,j-1], unchanged */
        float** bpp, /* neighbour entries for [i+1,j+1], unchanged */
        float** bmm, /* neighbour entries for [i-1,j-1], unchanged */
        float** bpm, /* neighbour entries for [i+1,j-1], unchanged */
        float** bmp, /* neighbour entries for [i-1,j+1], unchanged */
        float** f,   /* right hand side, unchanged */
        long* k,     /* number of iterations, output */
        float** u)   /* old and new solution, changed */

/*
   method of conjugate gradients without preconditioning for solving a
   linear system B u = f with a symmetric, nonadiagonal system matrix B
   that involves all eight 2D neighbours.
*/

{
    long i, j;   /* loop variables */
    float** r;   /* residue */
    float** p;   /* A-conjugate basis vector */
    float** q;   /* q = A p */
    float eps;   /* squared norm of residue r */
    float alpha; /* step size for updating u and r */
    float beta;  /* step size for updating p */
    float delta; /* auxiliary variable */
    float rho;   /* squared norm of right hand side f */

    /* ---- allocate storage ---- */

    alloc_matrix(&r, nx + 2, ny + 2);
    alloc_matrix(&p, nx + 2, ny + 2);
    alloc_matrix(&q, nx + 2, ny + 2);

    /* ---- INITIALISATIONS ---- */

    /* compute residue r = f - A * u */
    dummies_Neumann(u, nx, ny);
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
            r[i][j] = f[i][j]
                - (boo[i][j] * u[i][j] + bpo[i][j] * u[i + 1][j] + bmo[i][j] * u[i - 1][j]
                    + bop[i][j] * u[i][j + 1] + bom[i][j] * u[i][j - 1]
                    + bpp[i][j] * u[i + 1][j + 1] + bmm[i][j] * u[i - 1][j - 1]
                    + bpm[i][j] * u[i + 1][j - 1] + bmp[i][j] * u[i - 1][j + 1]);

    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
            p[i][j] = r[i][j];

    eps = inner_product(nx, ny, r, r);

    rho = inner_product(nx, ny, f, f);

    *k = 0;

    /* ---- ITERATIONS ---- */

    while ((eps >= 0.000000000000000001 * rho) && (*k < kmax))

    {

        /* compute q = A * p */
        matrix_times_vector(nx, ny, boo, bpo, bmo, bop, bom, bpp, bmm, bpm, bmp, p, q);

        /* update solution u and residue r */
        alpha = eps / inner_product(nx, ny, p, q);
        for (i = 1; i <= nx; i++)
            for (j = 1; j <= ny; j++) {
                u[i][j] = u[i][j] + alpha * p[i][j];
                r[i][j] = r[i][j] - alpha * q[i][j];
            }

        /* get next conjugate direction p */
        delta = eps;
        eps = inner_product(nx, ny, r, r);
        beta = eps / delta;
        for (i = 1; i <= nx; i++)
            for (j = 1; j <= ny; j++)
                p[i][j] = r[i][j] + beta * p[i][j];

        *k = *k + 1;

    } /* while */

    /* ---- disallocate storage ----*/

    disalloc_matrix(r, nx + 2, ny + 2);
    disalloc_matrix(p, nx + 2, ny + 2);
    disalloc_matrix(q, nx + 2, ny + 2);

    return;

} /* CG */

/*--------------------------------------------------------------------------*/

void GS

    (long kmax,      /* max. number of iterations */
        long nx,     /* image dimension in x direction */
        long ny,     /* image dimension in y direction */
        float** boo, /* diagonal entries for [i,j], unchanged */
        float** bpo, /* neighbour entries for [i+1,j], unchanged */
        float** bmo, /* neighbour entries for [i-1,j], unchanged */
        float** bop, /* neighbour entries for [i,j+1], unchanged */
        float** bom, /* neighbour entries for [i,j-1], unchanged */
        float** bpp, /* neighbour entries for [i+1,j+1], unchanged */
        float** bmm, /* neighbour entries for [i-1,j-1], unchanged */
        float** bpm, /* neighbour entries for [i+1,j-1], unchanged */
        float** bmp, /* neighbour entries for [i-1,j+1], unchanged */
        float** f,   /* right hand side, unchanged */
        long* k,     /* number of iterations, output */
        float** u)   /* old and new solution, changed */

/*
   Gauss-Seidel method for solving a linear system B u = f with a
   nonadiagonal system matrix B that involves all eight 2D neighbours.
*/

{
    long i, j; /* loop variables */
    float** r; /* residue */
    float eps; /* squared norm of residue r */
    float rho; /* squared norm of right hand side f */

    /* ---- allocate storage ---- */

    alloc_matrix(&r, nx + 2, ny + 2);

    /* ---- INITIALISATIONS ---- */

    /* ---- compute residue r = f - A * u ---- */

    dummies_Neumann(u, nx, ny);

    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
            r[i][j] = f[i][j]
                - (boo[i][j] * u[i][j] + bpo[i][j] * u[i + 1][j] + bmo[i][j] * u[i - 1][j]
                    + bop[i][j] * u[i][j + 1] + bom[i][j] * u[i][j - 1]
                    + bpp[i][j] * u[i + 1][j + 1] + bmm[i][j] * u[i - 1][j - 1]
                    + bpm[i][j] * u[i + 1][j - 1] + bmp[i][j] * u[i - 1][j + 1]);

    eps = inner_product(nx, ny, r, r);

    rho = inner_product(nx, ny, f, f);

    *k = 0;

    /* ---- ITERATIONS ---- */

    while ((eps >= 0.000000000000000001 * rho) && (*k < kmax))

    {
        /* ---- perform one Gauss-Seidel step ---- */

        dummies_Neumann(u, nx, ny);

        for (i = 1; i <= nx; i++)
            for (j = 1; j <= ny; j++)
                u[i][j] = (f[i][j] - bpo[i][j] * u[i + 1][j] - bmo[i][j] * u[i - 1][j]
                              - bop[i][j] * u[i][j + 1] - bom[i][j] * u[i][j - 1]
                              - bpp[i][j] * u[i + 1][j + 1] - bmm[i][j] * u[i - 1][j - 1]
                              - bpm[i][j] * u[i + 1][j - 1] - bmp[i][j] * u[i - 1][j + 1])
                    / boo[i][j];

        /* ---- compute residue r = f - A * u ---- */

        dummies_Neumann(u, nx, ny);

        for (i = 1; i <= nx; i++)
            for (j = 1; j <= ny; j++)
                r[i][j] = f[i][j]
                    - (boo[i][j] * u[i][j] + bpo[i][j] * u[i + 1][j] + bmo[i][j] * u[i - 1][j]
                        + bop[i][j] * u[i][j + 1] + bom[i][j] * u[i][j - 1]
                        + bpp[i][j] * u[i + 1][j + 1] + bmm[i][j] * u[i - 1][j - 1]
                        + bpm[i][j] * u[i + 1][j - 1] + bmp[i][j] * u[i - 1][j + 1]);

        eps = inner_product(nx, ny, r, r);

        /* ---- update k ---- */

        *k = *k + 1;

    } /* while */

    /* ---- disallocate storage ----*/

    disalloc_matrix(r, nx + 2, ny + 2);

    return;

} /* GS */

/*--------------------------------------------------------------------------*/

void eed_ex

    (float ht,        /* time step size */
        long nx,      /* image dimension in x direction */
        long ny,      /* image dimension in y direction */
        float hx,     /* pixel size in x direction */
        float hy,     /* pixel size in y direction */
        long dtype,   /* type of diffusivity */
        float lambda, /* contrast parameter */
        float sigma,  /* noise scale */
        float rho,    /* integration scale */
        float alpha,  /* dissipativity parameter */
        float gamma,  /* nonnegativity parameter */
        float** a,    /* confidence map, unchanged */
        float** u)    /* input: original image;  output: smoothed */

/*
 Edge-enhancing anisotropic diffusion filtering.
 2-stable explicit discretisation.
*/

{
    long i, j;                 /* loop variables */
    float** f;                 /* work copy of u */
    float **dxx, **dxy, **dyy; /* entries of structure/diffusion tensor */
    float** woo;               /* weights for [i,j] */
    float** wpo;               /* weights for [i+1,j] */
    float** wmo;               /* weights for [i-1,j] */
    float** wop;               /* weights for [i,j+1] */
    float** wom;               /* weights for [i,j-1] */
    float** wpp;               /* weights for [i+1,j+1] */
    float** wmm;               /* weights for [i-1,j-1] */
    float** wpm;               /* weights for [i+1,j-1] */
    float** wmp;               /* weights for [i-1,j+1] */

    /* ---- allocate storage ---- */

    alloc_matrix(&f, nx + 2, ny + 2);
    alloc_matrix(&dxx, nx + 1, ny + 1);
    alloc_matrix(&dxy, nx + 1, ny + 1);
    alloc_matrix(&dyy, nx + 1, ny + 1);
    alloc_matrix(&woo, nx + 2, ny + 2);
    alloc_matrix(&wpo, nx + 2, ny + 2);
    alloc_matrix(&wmo, nx + 2, ny + 2);
    alloc_matrix(&wop, nx + 2, ny + 2);
    alloc_matrix(&wom, nx + 2, ny + 2);
    alloc_matrix(&wpp, nx + 2, ny + 2);
    alloc_matrix(&wmm, nx + 2, ny + 2);
    alloc_matrix(&wpm, nx + 2, ny + 2);
    alloc_matrix(&wmp, nx + 2, ny + 2);

    /* ---- compute stencil weights ---- */

    /* copy u into f */
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
            f[i][j] = u[i][j];

    /* compute structure tensor on staggered grid (alters u!!!) */
    struct_tensor(u, nx, ny, hx, hy, sigma, rho, dxx, dxy, dyy);

    /* compute diffusion tensor on staggered grid */
    diff_tensor(dtype, lambda, nx - 1, ny - 1, dxx, dxy, dyy);

    /* compute stencil weights (on original grid) */
    weights(dxx, dxy, dyy, a, nx, ny, hx, hy, alpha, gamma, woo, wpo, wmo, wop, wom, wpp, wmm, wpm,
        wmp);

    /* ---- explicit diffusion ---- */

    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
            u[i][j] = f[i][j]
                + ht
                    * (woo[i][j] * f[i][j] + wpo[i][j] * f[i + 1][j] + wop[i][j] * f[i][j + 1]
                        + wpp[i][j] * f[i + 1][j + 1] + wpm[i][j] * f[i + 1][j - 1]
                        + wmo[i][j] * f[i - 1][j] + wom[i][j] * f[i][j - 1]
                        + wmm[i][j] * f[i - 1][j - 1] + wmp[i][j] * f[i - 1][j + 1]);

    /* ---- disallocate storage ---- */

    disalloc_matrix(f, nx + 2, ny + 2);
    disalloc_matrix(dxx, nx + 1, ny + 1);
    disalloc_matrix(dxy, nx + 1, ny + 1);
    disalloc_matrix(dyy, nx + 1, ny + 1);
    disalloc_matrix(woo, nx + 2, ny + 2);
    disalloc_matrix(wpo, nx + 2, ny + 2);
    disalloc_matrix(wmo, nx + 2, ny + 2);
    disalloc_matrix(wop, nx + 2, ny + 2);
    disalloc_matrix(wom, nx + 2, ny + 2);
    disalloc_matrix(wpp, nx + 2, ny + 2);
    disalloc_matrix(wmm, nx + 2, ny + 2);
    disalloc_matrix(wpm, nx + 2, ny + 2);
    disalloc_matrix(wmp, nx + 2, ny + 2);

    return;

} /* eed_ex */

/*--------------------------------------------------------------------------*/

void eed_im

    (float ht,        /* time step size */
        long nx,      /* image dimension in x direction */
        long ny,      /* image dimension in y direction */
        float hx,     /* pixel size in x direction */
        float hy,     /* pixel size in y direction */
        long dtype,   /* type of diffusivity */
        float lambda, /* contrast parameter */
        float sigma,  /* noise scale */
        float rho,    /* integration scale */
        float alpha,  /* dissipativity parameter */
        float gamma,  /* nonnegativity parameter */
        float stype,  /* type of linear solver */
        long imax,    /* max. number of solver iterations */
        float** a,    /* confidence map, unchanged */
        long* count,  /* number of linear solver iterations, output */
        float** u)    /* input: original image;  output: smoothed */

/*
 Edge-enhancing anisotropic diffusion filtering.
 2-stable semi-implicit discretisation.
*/

{
    long i, j;                 /* loop variables */
    float** f;                 /* work copy of u */
    float** v;                 /* intermediate result */
    float **dxx, **dxy, **dyy; /* entries of structure/diffusion tensor */
    float** woo;               /* weights for [i,j] */
    float** wpo;               /* weights for [i+1,j] */
    float** wmo;               /* weights for [i-1,j] */
    float** wop;               /* weights for [i,j+1] */
    float** wom;               /* weights for [i,j-1] */
    float** wpp;               /* weights for [i+1,j+1] */
    float** wmm;               /* weights for [i-1,j-1] */
    float** wpm;               /* weights for [i+1,j-1] */
    float** wmp;               /* weights for [i-1,j+1] */

    /* ---- allocate storage ---- */

    alloc_matrix(&f, nx + 2, ny + 2);
    alloc_matrix(&v, nx + 2, ny + 2);
    alloc_matrix(&dxx, nx + 1, ny + 1);
    alloc_matrix(&dxy, nx + 1, ny + 1);
    alloc_matrix(&dyy, nx + 1, ny + 1);
    alloc_matrix(&woo, nx + 2, ny + 2);
    alloc_matrix(&wpo, nx + 2, ny + 2);
    alloc_matrix(&wmo, nx + 2, ny + 2);
    alloc_matrix(&wop, nx + 2, ny + 2);
    alloc_matrix(&wom, nx + 2, ny + 2);
    alloc_matrix(&wpp, nx + 2, ny + 2);
    alloc_matrix(&wmm, nx + 2, ny + 2);
    alloc_matrix(&wpm, nx + 2, ny + 2);
    alloc_matrix(&wmp, nx + 2, ny + 2);

    /* ---- compute stencil weights ---- */

    /* copy u into f */
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
            f[i][j] = u[i][j];

    /* compute structure tensor on staggered grid (alters u!!!) */
    struct_tensor(u, nx, ny, hx, hy, sigma, rho, dxx, dxy, dyy);

    /* compute diffusion tensor on staggered grid */
    diff_tensor(dtype, lambda, nx - 1, ny - 1, dxx, dxy, dyy);

    /* compute stencil weights (on original grid) */
    weights(dxx, dxy, dyy, a, nx, ny, hx, hy, alpha, gamma, woo, wpo, wmo, wop, wom, wpp, wmm, wpm,
        wmp);

    /* ---- step 1:  u := A * f ---- */

    dummies_Neumann(f, nx, ny);
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
            u[i][j] = woo[i][j] * f[i][j] + wpo[i][j] * f[i + 1][j] + wmo[i][j] * f[i - 1][j]
                + wop[i][j] * f[i][j + 1] + wom[i][j] * f[i][j - 1] + wpp[i][j] * f[i + 1][j + 1]
                + wmm[i][j] * f[i - 1][j - 1] + wpm[i][j] * f[i + 1][j - 1]
                + wmp[i][j] * f[i - 1][j + 1];

    /* ---- step 2:  solve (I - ht * A) v = u for v ---- */

    /* compute system matrix (I - ht * A) */
    weights_to_matrix(nx, ny, ht, woo, wpo, wmo, wop, wom, wpp, wmm, wpm, wmp);

    /* initialise v */
    for (i = 0; i <= nx + 1; i++)
        for (j = 0; j <= ny + 1; j++)
            v[i][j] = 0.0;

    /* iterative solver with at most 10000 iterations */
    if (stype == 0)
        /* conjugate gradients */
        CG(imax, nx, ny, woo, wpo, wmo, wop, wom, wpp, wmm, wpm, wmp, u, &*count, v);
    else if (stype == 1)
        /* Gauss-Seidel */
        GS(imax, nx, ny, woo, wpo, wmo, wop, wom, wpp, wmm, wpm, wmp, u, &*count, v);

    /* ---- step 3:  u = f + ht * v ---- */

    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
            u[i][j] = f[i][j] + ht * v[i][j];

    /* ---- disallocate storage ---- */

    disalloc_matrix(f, nx + 2, ny + 2);
    disalloc_matrix(v, nx + 2, ny + 2);
    disalloc_matrix(dxx, nx + 1, ny + 1);
    disalloc_matrix(dxy, nx + 1, ny + 1);
    disalloc_matrix(dyy, nx + 1, ny + 1);
    disalloc_matrix(woo, nx + 2, ny + 2);
    disalloc_matrix(wpo, nx + 2, ny + 2);
    disalloc_matrix(wmo, nx + 2, ny + 2);
    disalloc_matrix(wop, nx + 2, ny + 2);
    disalloc_matrix(wom, nx + 2, ny + 2);
    disalloc_matrix(wpp, nx + 2, ny + 2);
    disalloc_matrix(wmm, nx + 2, ny + 2);
    disalloc_matrix(wpm, nx + 2, ny + 2);
    disalloc_matrix(wmp, nx + 2, ny + 2);

    return;

} /* eed_im */

/*--------------------------------------------------------------------------*/

float elliptic_residue

    (long nx,         /* image dimension in x direction */
        long ny,      /* image dimension in y direction */
        float hx,     /* pixel size in x direction */
        float hy,     /* pixel size in y direction */
        long dtype,   /* type of diffusivity */
        float lambda, /* contrast parameter */
        float sigma,  /* noise scale */
        float rho,    /* integration scale */
        float alpha,  /* dissipativity parameter */
        float gamma,  /* nonnegativity parameter */
        float** a,    /* confidence map, unchanged */
        float** u)    /* solution, unchanged */

/*
 computes the normalised 2-norm of the residue of the elliptic equation
*/

{
    long i, j;                 /* loop variables */
    float res;                 /* normalised 2-norm of elliptic residue */
    float** f;                 /* work copy of u */
    float** b;                 /* right hand side */
    float **dxx, **dxy, **dyy; /* entries of structure/diffusion tensor */
    float** woo;               /* weights for [i,j] */
    float** wpo;               /* weights for [i+1,j] */
    float** wmo;               /* weights for [i-1,j] */
    float** wop;               /* weights for [i,j+1] */
    float** wom;               /* weights for [i,j-1] */
    float** wpp;               /* weights for [i+1,j+1] */
    float** wmm;               /* weights for [i-1,j-1] */
    float** wpm;               /* weights for [i+1,j-1] */
    float** wmp;               /* weights for [i-1,j+1] */

    /* ---- allocate storage ---- */

    alloc_matrix(&f, nx + 2, ny + 2);
    alloc_matrix(&b, nx + 2, ny + 2);
    alloc_matrix(&dxx, nx + 1, ny + 1);
    alloc_matrix(&dxy, nx + 1, ny + 1);
    alloc_matrix(&dyy, nx + 1, ny + 1);
    alloc_matrix(&woo, nx + 2, ny + 2);
    alloc_matrix(&wpo, nx + 2, ny + 2);
    alloc_matrix(&wmo, nx + 2, ny + 2);
    alloc_matrix(&wop, nx + 2, ny + 2);
    alloc_matrix(&wom, nx + 2, ny + 2);
    alloc_matrix(&wpp, nx + 2, ny + 2);
    alloc_matrix(&wmm, nx + 2, ny + 2);
    alloc_matrix(&wpm, nx + 2, ny + 2);
    alloc_matrix(&wmp, nx + 2, ny + 2);

    /* ---- compute stencil weights ---- */

    /* copy u into f */
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
            f[i][j] = u[i][j];

    /* compute structure tensor on staggered grid (alters f!!!) */
    struct_tensor(f, nx, ny, hx, hy, sigma, rho, dxx, dxy, dyy);

    /* compute diffusion tensor on staggered grid */
    diff_tensor(dtype, lambda, nx - 1, ny - 1, dxx, dxy, dyy);

    /* compute stencil weights (on original grid) */
    weights(dxx, dxy, dyy, a, nx, ny, hx, hy, alpha, gamma, woo, wpo, wmo, wop, wom, wpp, wmm, wpm,
        wmp);

    /* modify stencil weights and compute right hand side */
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
            if (a[i][j] > 0.0) {
                woo[i][j] = 1.0;
                b[i][j] = u[i][j];
            } else
                b[i][j] = 0.0;

    /* ---- compute residue ---- */

    res = residue(nx, ny, woo, wpo, wmo, wop, wom, wpp, wmm, wpm, wmp, b, u);

    /* ---- disallocate storage ---- */

    disalloc_matrix(f, nx + 2, ny + 2);
    disalloc_matrix(b, nx + 2, ny + 2);
    disalloc_matrix(dxx, nx + 1, ny + 1);
    disalloc_matrix(dxy, nx + 1, ny + 1);
    disalloc_matrix(dyy, nx + 1, ny + 1);
    disalloc_matrix(woo, nx + 2, ny + 2);
    disalloc_matrix(wpo, nx + 2, ny + 2);
    disalloc_matrix(wmo, nx + 2, ny + 2);
    disalloc_matrix(wop, nx + 2, ny + 2);
    disalloc_matrix(wom, nx + 2, ny + 2);
    disalloc_matrix(wpp, nx + 2, ny + 2);
    disalloc_matrix(wmm, nx + 2, ny + 2);
    disalloc_matrix(wpm, nx + 2, ny + 2);
    disalloc_matrix(wmp, nx + 2, ny + 2);

    return (res);

} /* elliptic_residue */

/*--------------------------------------------------------------------------*/

void analyse

    (float** u,       /* image, unchanged */
        long nx,      /* pixel number in x direction */
        long ny,      /* pixel number in x direction */
        float* min,   /* minimum, output */
        float* max,   /* maximum, output */
        float* mean,  /* mean, output */
        float* stdev, /* standard deviation, output */
        float* norm)  /* L2 norm, output */

/*
 calculates minimum, maximum, mean, standard deviation, and norm of an
 image u
*/

{
    long i, j;    /* loop variables */
    float help;   /* auxiliary variable */
    double help2; /* auxiliary variable */

    *min = u[1][1];
    *max = u[1][1];
    help2 = 0.0;
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++) {
            if (u[i][j] < *min)
                *min = u[i][j];
            if (u[i][j] > *max)
                *max = u[i][j];
            help2 = help2 + (double)u[i][j];
        }
    *mean = (float)help2 / (nx * ny);

    *norm = 0.0;
    *stdev = 0.0;
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++) {
            *norm = *norm + u[i][j] * u[i][j];
            help = u[i][j] - *mean;
            *stdev = *stdev + help * help;
        }
    *stdev = sqrt(*stdev / (nx * ny));
    *norm = sqrt(*norm / (nx * ny));
    return;

} /* analyse */

/*--------------------------------------------------------------------------*/

void inpainting(char* in, char* mask, char* out, long dtype, float lambda, float sigma, float rho,
    float alpha, float gamma, long timediscr, float ht, long kmax, long stype, long imax)

{
    char row[80];             /* for reading data */
    float** f;                /* evolving image */
    float** a;                /* inpainting mask, 0 for missing data */
    long i, j, k;             /* loop variables */
    long nx, ny;              /* image size in x, y direction */
    FILE *inimage, *outimage; /* input file, output file */
    long count;               /* number of linear solver iterations */
    float hx;                 /* step size in x direction */
    float hy;                 /* step size in y direction */
    float max, min;           /* largest, smallest grey value */
    float mean;               /* average grey value */
    float stdev;              /* standard deviation */
    float norm;               /* normalised 2-norm of the solution */
    float res;                /* normalised 2-norm of elliptic residue */
    unsigned char byte;       /* for data conversion */

    printf("Parameter values:\n\tdtype: %ld\n\tlambda: %lf\n\tsigma: %lf\n\trho: %lf\n\talpha: "
           "%lf\n\tgamma: %lf\n\tht: %lf\n\tkmax: %ld\n\tstype: %ld\n\timax: %ld",
        dtype, lambda, sigma, rho, alpha, gamma, ht, kmax, stype, imax);

    printf("\n");
    printf("EED-BASED INPAINTING\n\n");
    printf("***************************************************\n\n");
    printf("    Copyright 2008 by Joachim Weickert             \n");
    printf("    Faculty of Mathematics and Computer Science    \n");
    printf("    Saarland University, Germany                   \n\n");
    printf("    All rights reserved. Unauthorized usage,       \n");
    printf("    copying, hiring, and selling prohibited.       \n\n");
    printf("    Send bug reports to                            \n");
    printf("    weickert@mia.uni-saarland.de                   \n\n");
    printf("***************************************************\n\n");

    /* open pgm file and read header */
    inimage = fopen(in, "r");
    fgets(row, 300, inimage);
    fgets(row, 300, inimage);
    while (row[0] == '#')
        fgets(row, 300, inimage);
    sscanf(row, "%ld %ld", &nx, &ny);
    fgets(row, 300, inimage);

    /* allocate storage */
    alloc_matrix(&f, nx + 2, ny + 2);

    /* read image data */
    for (j = 1; j <= ny; j++)
        for (i = 1; i <= nx; i++)
            f[i][j] = (float)getc(inimage);
    fclose(inimage);

    /* ---- read inpainting mask (pgm format P5) ---- */

    /* open pgm file and read header */
    inimage = fopen(mask, "r");
    fgets(row, 300, inimage);
    fgets(row, 300, inimage);
    while (row[0] == '#')
        fgets(row, 300, inimage);
    sscanf(row, "%ld %ld", &nx, &ny);
    fgets(row, 300, inimage);

    /* allocate storage */
    alloc_matrix(&a, nx + 2, ny + 2);

    /* read image data */
    for (j = 1; j <= ny; j++)
        for (i = 1; i <= nx; i++)
            a[i][j] = (float)getc(inimage);
    fclose(inimage);

    /* ---- initialisations ---- */

    /* normalise inpainting mask */
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
            a[i][j] = a[i][j] / 255.0;

    /* unit grid size */
    hx = hy = 1.0;

    /* initialise loop counter */
    k = 0;

    /* check minimum, maximum, mean, variance and residue */

    analyse(f, nx, ny, &min, &max, &mean, &stdev, &norm);
    res = elliptic_residue(nx, ny, hx, hy, dtype, lambda, sigma, rho, alpha, gamma, a, f);

    printf("initial image\n");
    printf("minimum:         %8.2lf \n", min);
    printf("maximum:         %8.2lf \n", max);
    printf("mean:            %8.2lf \n", mean);
    printf("std. dev.:       %8.2lf \n", stdev);
    printf("2-norm:          %8.2lf \n", norm);
    printf("residue:         %8.6lf \n\n", res);

    /* ---- process image ---- */

    while (k < kmax) {
        /* perform one inpainting iteration */
        k = k + 1;
        if (timediscr == 0)
            /* explicit scheme */
            eed_ex(ht, nx, ny, hx, hy, dtype, lambda, sigma, rho, alpha, gamma, a, f);
        else if (timediscr == 1)
            /* semi-implicit scheme */
            eed_im(ht, nx, ny, hx, hy, dtype, lambda, sigma, rho, alpha, gamma, stype, imax, a,
                &count, f);

        /* check stable pixels and other things every 10th iteration */
        if (10 * (k / 10) == k) {
            analyse(f, nx, ny, &min, &max, &mean, &stdev, &norm);
            res = elliptic_residue(nx, ny, hx, hy, dtype, lambda, sigma, rho, alpha, gamma, a, f);

            printf("iteration number: %7ld \n", k);
            if (timediscr == 1)
                printf("inner iterations: %7ld \n", count);
            printf("minimum:         %8.2lf \n", min);
            printf("maximum:         %8.2lf \n", max);
            printf("mean:            %8.2lf \n", mean);
            printf("std. dev.:       %8.2lf \n", stdev);
            printf("2-norm:          %8.2lf \n", norm);
            printf("residue:         %8.6lf \n\n", res);
        }
    }

    /* ---- write output image (pgm format P5) ---- */

    /* open file and write header (incl. filter parameters) */
    outimage = fopen(out, "w");
    fprintf(outimage, "P5 \n");
    fprintf(outimage, "# inpainting with EED\n");
    fprintf(outimage, "# semi-implicit scheme\n");
    if (timediscr == 0)
        fprintf(outimage, "# explicit scheme\n");
    else if (timediscr == 1) {
        fprintf(outimage, "# semi-implicit scheme\n");
        if (stype == 0)
            fprintf(outimage, "# linear solver:    conjugate gradients\n");
        else if (stype == 1)
            fprintf(outimage, "# linear solver:    Gauss-Seidel\n");
        fprintf(outimage, "# imax:             %1ld\n", imax);
    }
    fprintf(outimage, "# initial image:    %s\n", in);
    fprintf(outimage, "# inpainting mask:  %s\n", mask);
    if (dtype == 0)
        fprintf(outimage, "# diffusivity:      Charbonnier\n");
    else if (dtype == 1)
        fprintf(outimage, "# diffusivity:      Weickert\n");
    fprintf(outimage, "# lambda:           %8.6lf\n", lambda);
    fprintf(outimage, "# sigma:          %8.4lf\n", sigma);
    fprintf(outimage, "# rho:            %8.4lf\n", rho);
    fprintf(outimage, "# alpha:          %8.4lf\n", alpha);
    fprintf(outimage, "# gamma:          %8.4lf\n", gamma);
    fprintf(outimage, "# ht:             %8.2lf\n", ht);
    fprintf(outimage, "# time steps:     %8ld\n", kmax);
    fprintf(outimage, "# minimum:        %8.2lf\n", min);
    fprintf(outimage, "# maximum:        %8.2lf\n", max);
    fprintf(outimage, "# mean:           %8.2lf\n", mean);
    fprintf(outimage, "# std. dev.:      %8.2lf\n", stdev);
    fprintf(outimage, "# 2-norm:         %8.2lf\n", norm);
    fprintf(outimage, "# residue:        %8.6lf\n", res);
    fprintf(outimage, "%ld %ld \n255\n", nx, ny);

    /* write image data and close file */
    for (j = 1; j <= ny; j++)
        for (i = 1; i <= nx; i++) {
            if (f[i][j] < 0.0)
                byte = (unsigned char)(0.0);
            else if (f[i][j] > 255.0)
                byte = (unsigned char)(255.0);
            else
                byte = (unsigned char)(f[i][j] + 0.499999);
            fwrite(&byte, sizeof(unsigned char), 1, outimage);
        }
    fclose(outimage);
    printf("output image %s successfully written\n\n", out);

    /* ---- disallocate storage ---- */

    disalloc_matrix(f, nx + 2, ny + 2);
    disalloc_matrix(a, nx + 2, ny + 2);
}

/*--------------------------------------------------------------------------*/

/*
    inpainting with default values as described in paramerers
*/
void inpainting_compression(char* in, char* mask, char* out, float lambda, float sigma,
    float alpha, float gamma, long outer_iter, long solver_iter)
{
    inpainting(
        in, mask, out, 0, lambda, sigma, 0.0, alpha, gamma, 1, 1000, outer_iter, 0, solver_iter);
}

/*--------------------------------------------------------------------------*/

float MSE(float** u, float** v, long nx, long ny)
{
    float sum = 0;
    for (long i = 1; i <= nx; ++i) {
        for (long j = 1; j <= ny; ++j) {
            sum += pow(u[i][j] - v[i][j], 2);
        }
    }
    return sum / (nx * ny);
}

float PSNR(float** u, float** v, long nx, long ny)
{
    float max = -INFINITY;
    float mse = MSE(u, v, nx, ny);
    for (long i = 1; i <= nx; ++i) {
        for (long j = 1; j <= ny; ++j) {
            if (u[i][j] > max)
                max = u[i][j];
        }
    }
    return 20 * log10(max) - 10 * log10(mse);
}

/*--------------------------------------------------------------------------*/
int main(int argc, char** argv)

{
    char row[80];             /* for reading data */
    char* in;                 /* for reading data */
    char* in2;                /* for reading data */
    float** f;                /* evolving image */
    float** u;                /* Copy of initial image */
    float** a;                /* inpainting mask, 0 for missing data */
    long i, j, k;             /* loop variables */
    long nx, ny;              /* image size in x, y direction */
    FILE *inimage, *outimage; /* input file, output file */
    long count;               /* number of linear solver iterations */
    float hx;                 /* step size in x direction */
    float hy;                 /* step size in y direction */
    float max, min;           /* largest, smallest grey value */
    float mean;               /* average grey value */
    float stdev;              /* standard deviation */
    float norm;               /* normalised 2-norm of the solution */
    float res;                /* normalised 2-norm of elliptic residue */
    unsigned char byte;       /* for data conversion */

    /* Defaults or command line arguments */
    char* out = "inpaint.pgm";                /* for reading data */
    long dtype = 0;      /* type of diffusivity */
    float lambda = 0.03; /* contrast parameter */
    float sigma = 3.0;   /* noise scale */
    float rho = 0;       /* integration scale */
    float alpha = 0.5;   /* dissipativity parameter */
    float gamma = 1;     /* nonnegativity parameter */
    float ht = 1000;     /* time step size */
    long kmax = 100;     /* largest iteration number */
    long imax = 200;     /* largest solver iteration number */
    long timediscr = 1;  /* type of time discretisation */
    long stype = 0;      /* type of linear system solver */

    int c;
    while ((c = getopt(argc, argv, "l:s:a:g:n:N:o:")) != -1) {
        switch (c) {
        case 'l':
            lambda = atof(optarg);
            break;
        case 's':
            sigma = atof(optarg);
            break;
        case 'a':
            alpha = atof(optarg);
            break;
        case 'g':
            gamma = atof(optarg);
            break;
        case 'n':
            kmax = atol(optarg);
            break;
        case 'N':
            imax = atol(optarg);
            break;
        case 'o':
            out = optarg;
            break;
        default:
            printf("Error reading arguments!\n");
            abort();
        }
    }

    if (optind == argc) {
        printf("No input files specified. Aborting...\n");
        abort();
    }

    if (argc - optind > 2) {
        printf("Too many positional arguments. Aborting...\n");
        abort();
    }

    if (argc - optind < 2) {
        printf("No mask specified. Aborting...\n");
        abort();
    }

    in = argv[optind];
    in2 = argv[optind + 1];

    /* ---- read input image (pgm format P5) ---- */

    /* open pgm file and read header */
    inimage = fopen(in, "r");
    fgets(row, 300, inimage);
    fgets(row, 300, inimage);
    while (row[0] == '#')
        fgets(row, 300, inimage);
    sscanf(row, "%ld %ld", &nx, &ny);
    fgets(row, 300, inimage);

    /* allocate storage */
    alloc_matrix(&f, nx + 2, ny + 2);

    /* read image data */
    for (j = 1; j <= ny; j++)
        for (i = 1; i <= nx; i++)
            f[i][j] = (float)getc(inimage);
    fclose(inimage);

    /* Copy f into u for MSE calculation */
    alloc_matrix(&u, nx + 2, ny + 2);
    for (j = 1; j <= ny; j++)
        for (i = 1; i <= nx; i++)
            u[i][j] = f[i][j];

    /* ---- read inpainting mask (pgm format P5) ---- */

    /* open pgm file and read header */
    inimage = fopen(in2, "r");
    fgets(row, 300, inimage);
    fgets(row, 300, inimage);
    while (row[0] == '#')
        fgets(row, 300, inimage);
    sscanf(row, "%ld %ld", &nx, &ny);
    fgets(row, 300, inimage);

    /* allocate storage */
    alloc_matrix(&a, nx + 2, ny + 2);

    /* read image data */
    for (j = 1; j <= ny; j++)
        for (i = 1; i <= nx; i++)
            a[i][j] = (float)getc(inimage);
    fclose(inimage);
    /* ---- initialisations ---- */

    /* normalise inpainting mask */
    for (i = 1; i <= nx; i++)
        for (j = 1; j <= ny; j++)
            a[i][j] = a[i][j] / 255.0;

    /* unit grid size */
    hx = hy = 1.0;

    /* initialise loop counter */
    k = 0;

    printf("\n");
    printf("EED-BASED INPAINTING\n\n");
    printf("***************************************************\n\n");
    printf("    Copyright 2008 by Joachim Weickert             \n");
    printf("    Faculty of Mathematics and Computer Science    \n");
    printf("    Saarland University, Germany                   \n\n");
    printf("    All rights reserved. Unauthorized usage,       \n");
    printf("    copying, hiring, and selling prohibited.       \n\n");
    printf("    Send bug reports to                            \n");
    printf("    weickert@mia.uni-saarland.de                   \n\n");
    printf("***************************************************\n\n");

    /* check minimum, maximum, mean, variance and residue */

    analyse(f, nx, ny, &min, &max, &mean, &stdev, &norm);
    res = elliptic_residue(nx, ny, hx, hy, dtype, lambda, sigma, rho, alpha, gamma, a, f);

    printf("initial image\n");
    printf("minimum:         %8.2lf \n", min);
    printf("maximum:         %8.2lf \n", max);
    printf("mean:            %8.2lf \n", mean);
    printf("std. dev.:       %8.2lf \n", stdev);
    printf("2-norm:          %8.2lf \n", norm);
    printf("residue:         %8.6lf \n\n", res);

    /* ---- process image ---- */

    while (res >= 0.000001 && k < kmax) {
        /* perform one inpainting iteration */
        k = k + 1;
        if (timediscr == 0)
            /* explicit scheme */
            eed_ex(ht, nx, ny, hx, hy, dtype, lambda, sigma, rho, alpha, gamma, a, f);
        else if (timediscr == 1)
            /* semi-implicit scheme */
            eed_im(ht, nx, ny, hx, hy, dtype, lambda, sigma, rho, alpha, gamma, stype, imax, a,
                &count, f);

        /* check stable pixels and other things every 10th iteration */
        if (10 * (k / 10) == k) {
            analyse(f, nx, ny, &min, &max, &mean, &stdev, &norm);
            res = elliptic_residue(nx, ny, hx, hy, dtype, lambda, sigma, rho, alpha, gamma, a, f);

            printf("iteration number: %7ld \n", k);
            if (timediscr == 1)
                printf("inner iterations: %7ld \n", count);
            printf("minimum:         %8.2lf \n", min);
            printf("maximum:         %8.2lf \n", max);
            printf("mean:            %8.2lf \n", mean);
            printf("std. dev.:       %8.2lf \n", stdev);
            printf("2-norm:          %8.2lf \n", norm);
            printf("residue:         %8.6lf \n\n", res);
        } /* if */
    }     /* while */

    float mse = MSE(u, f, nx, ny);
    printf("MSE between inpaint and original image: %f\n", mse);

    float psnr = PSNR(u, f, nx, ny);
    printf("PSNR between inpaint and original image: %f\n", psnr);

    /* ---- write output image (pgm format P5) ---- */

    /* open file and write header (incl. filter parameters) */
    outimage = fopen(out, "w");
    fprintf(outimage, "P5 \n");
    fprintf(outimage, "# inpainting with EED\n");
    if (timediscr == 0)
        fprintf(outimage, "# explicit scheme\n");
    else if (timediscr == 1) {
        fprintf(outimage, "# semi-implicit scheme\n");
        if (stype == 0)
            fprintf(outimage, "# linear solver:    conjugate gradients\n");
        else if (stype == 1)
            fprintf(outimage, "# linear solver:    Gauss-Seidel\n");
        fprintf(outimage, "# imax:             %1ld\n", imax);
    }
    fprintf(outimage, "# initial image:    %s\n", in);
    fprintf(outimage, "# inpainting mask:  %s\n", in2);
    if (dtype == 0)
        fprintf(outimage, "# diffusivity:      Charbonnier\n");
    else if (dtype == 1)
        fprintf(outimage, "# diffusivity:      Weickert\n");
    fprintf(outimage, "# lambda:           %8.6lf\n", lambda);
    fprintf(outimage, "# sigma:          %8.4lf\n", sigma);
    fprintf(outimage, "# rho:            %8.4lf\n", rho);
    fprintf(outimage, "# alpha:          %8.4lf\n", alpha);
    fprintf(outimage, "# gamma:          %8.4lf\n", gamma);
    fprintf(outimage, "# ht:             %8.2lf\n", ht);
    fprintf(outimage, "# time steps:     %8ld\n", kmax);
    fprintf(outimage, "# minimum:        %8.2lf\n", min);
    fprintf(outimage, "# maximum:        %8.2lf\n", max);
    fprintf(outimage, "# mean:           %8.2lf\n", mean);
    fprintf(outimage, "# std. dev.:      %8.2lf\n", stdev);
    fprintf(outimage, "# 2-norm:         %8.2lf\n", norm);
    fprintf(outimage, "# residue:        %8.6lf\n", res);
    fprintf(outimage, "# MSE:            %8.6lf\n", mse);
    fprintf(outimage, "# PSNR:           %8.6lf\n", psnr);
    fprintf(outimage, "%ld %ld \n255\n", nx, ny);

    /* write image data and close file */
    for (j = 1; j <= ny; j++)
        for (i = 1; i <= nx; i++) {
            if (f[i][j] < 0.0)
                byte = (unsigned char)(0.0);
            else if (f[i][j] > 255.0)
                byte = (unsigned char)(255.0);
            else
                byte = (unsigned char)(f[i][j] + 0.499999);
            fwrite(&byte, sizeof(unsigned char), 1, outimage);
        }
    fclose(outimage);
    printf("output image %s successfully written\n\n", out);

    /* ---- disallocate storage ---- */

    disalloc_matrix(f, nx + 2, ny + 2);
    disalloc_matrix(a, nx + 2, ny + 2);
    return (0);
}
