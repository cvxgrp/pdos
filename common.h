#ifndef COMMON_H_GUARD
#define COMMON_H_GUARD

#include "pdos.h"
#include "cs.h"
#include <math.h>

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#define RATIO 1e6
#define SQRT_RATIO 1e3

/*
 * commonWorkInit(const Data *d, const Cone *k)
 * --------------------------------------------
 * Given a piece of data and a cone, this function performs the workspace
 * initalization that is *common* to both the direct and indirect methods.
 *
 * It will allocate memory for the PDOS algorithm: the vectors, x, stilde, s,
 * and y. Furthermore, x and stilde are contiguous in memory. It will allocate
 * memory for the scaling matrices D and E. It will also set the lambda
 * parameter.
 *
 * If we are to normalize the matrix A, a normalized copy of the data is made.
 *
 */

static inline Work *commonWorkInit(const Data *d, const Cone *k) {
  idxint i,j = 0;
  idxint Anz = d->Ap[d->n];

  Work * w = PDOS_malloc(sizeof(Work));
  // copy dimensions
  w->m = d->m; w->n = d->n;
  // copy parameters (pointer)
  w->params = d->p;
  // ensure that x, stilde are contiguous in memory
  w->x = PDOS_calloc(d->n + MAX(d->m,d->n),sizeof(double));
  if (d->x != NULL)
  {
    memcpy(w->x, d->x, d->n*sizeof(double));
  }
  w->stilde = w->x + d->n;
  // allocate workspace memory for s and y
  w->s = PDOS_calloc(d->m,sizeof(double));
  if (d->s != NULL)
  {
    memcpy(w->s, d->s, d->m*sizeof(double));
  }
  w->y = PDOS_calloc(d->m,sizeof(double));
  if (d->y != NULL)
  {
    memcpy(w->y, d->y, d->m*sizeof(double));
  }

  // allocate workspace memory for normalization matrices
  w->D = PDOS_calloc(d->m, sizeof(double));
  w->E = PDOS_calloc(d->n, sizeof(double));

  if(d->p->NORMALIZE) {
    w->Ax = PDOS_calloc(Anz, sizeof(double));
    w->Ai = d->Ai;
    w->Ap = d->Ap;
    w->b = PDOS_calloc(d->m, sizeof(double));
    w->c = PDOS_calloc(d->n, sizeof(double));

    // set w->Ax = d->Ax
    memcpy(w->Ax, d->Ax, Anz*sizeof(double));

    // set w->c = d->c
    memcpy(w->c, d->c, d->n*sizeof(double));

    // set w->b = d->b
    memcpy(w->b, d->b, d->m*sizeof(double));

    idxint ind, cone;
    // compute norm across rows
    for(i = 0; i < d->n; ++i) { // cols
      for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
        w->D[d->Ai[j]] += d->Ax[j] * d->Ax[j];
      }
    }
    // now collapse cones together by taking average norm square
    ind = k->f + k->l;
    double sum = 0;
    for(cone = 0; cone < k->qsize; ++cone) {
      sum = 0;
      for(i = 0; i < k->q[cone]; ++i) {
        sum += w->D[ind + i];
      }
      for(i = 0; i < k->q[cone]; ++i) {
        w->D[ind+i] = sum / k->q[cone];
      }
      ind += k->q[cone];
    }
    
    // walk through d->m rows of w->D, normalize and store in D
    for(i = 0; i < d->m; ++i) {
      w->D[i] = fabs(w->D[i]) > 1e-6 ? 1.0 / sqrt(w->D[i]) : 1.0;
    }
    
    // now scale A
    for(i = 0; i < d->n; ++i) { // cols
      for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
        w->Ax[j] *= w->D[d->Ai[j]];
      }
    }

    // now compute norms of columns (of normalized data)
    for(i = 0; i < d->n; ++i) { // cols
      for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
        w->E[i] += w->Ax[j] * w->Ax[j];
      }
    }

    for(i = 0; i < d->n; ++i) {
      w->E[i] = fabs(w->E[i]) > 1e-6 ? 1.0 / sqrt(w->E[i]) : 1.0;
    }

    // now scale A again
    for(i = 0; i < d->n; ++i) { // cols
      for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
        w->Ax[j] *= w->E[i];
      }
    }

    // scale c
    for(i = 0; i < d->n; ++i) {
      w->c[i] = d->c[i]*w->E[i];
    }
    // scale b
    for(i = 0; i < d->m; ++i) {
      w->b[i] = d->b[i]*w->D[i];
    }
    
    // y = D^{-1}*y (scale the initial variable)
    for(i = 0; i < w->m; ++i) {
      w->y[i] /= w->D[i];
    }

    // x = E^{-1}*x (scale the initial variable)
    for(i = 0; i < w->n; ++i) {
      w->x[i] /= w->E[i];      
    }

    // s = D*s (scale the initial variable)
    for(i = 0; i < w->m; ++i) {
      w->s[i] *= w->D[i];
    }

  } else {
    // if we don't normalize, we just point out workspace copy to the actual
    // data copies of the problem data
    w->Ax = d->Ax;
    w->Ai = d->Ai;
    w->Ap = d->Ap;
    w->b = d->b;
    w->c = d->c;

    idxint i;

    // set the scaling matrices to the identity
    for( i=0; i < d->m; ++i ) w->D[i] = 1.0;
    for( i=0; i < d->n; ++i ) w->E[i] = 1.0;
  }

  // PDOS_printf("||b||_2: %f ||c||_2: %f\n", calcNorm(w->b,w->m), calcNorm(w->c,w->n));
  w->lambda = (1e-6 + calcNorm(w->b,w->m)) / (1e-6 + calcNorm(w->c,w->n)) ;

  // set ratio of "x" space penalty (1e-6) to "s,y" space penalty (1)
  for( i=0; i < d->n; ++i ) {
    w->E[i] *= SQRT_RATIO;
    w->c[i] *= SQRT_RATIO;
    w->x[i] /= SQRT_RATIO;
  }
  for( i=0; i < Anz; ++i ) {
    w->Ax[i] *= SQRT_RATIO;
  }
  
  // transpose the A matrix and store it
  // first, store "A" in "cs" format
  cs * A = PDOS_calloc(1, sizeof(cs));
  A->m = w->m ;
  A->n = w->n ;
  A->nzmax = MAX (Anz, 1) ;
  A->nz = -1 ; // in compressed column form
  A->p = w->Ap;
  A->i = w->Ai;
  A->x = w->Ax;

  // now transpose
  cs * At = cs_transpose(A, 1);

  w->Atx = At->x;
  w->Ati = At->i;
  w->Atp = At->p;
  PDOS_free(At);  // "orphans" At->x, At->i, At->p, but they have been handed
                  // off to w->Atx, w->Ati, and w->Atp
  PDOS_free(A);

  return w;
}

#endif
