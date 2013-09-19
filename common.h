#ifndef COMMON_H_GUARD
#define COMMON_H_GUARD

#include "pdos.h"
#include "cs.h"
#include <math.h>

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#define EQUILIBRATE_ITERS 3
#define RATIO 1e6
#define SQRT_RATIO 1e3

static inline void collapseWorkspaceData(Work *w, const Data *d, const Cone *k)
{
  // thie function replaces all nonzero entires of w->Ax with max of absolute
  // value of each cone, same with w->b
  //
  // assumes w->Ax contains all of d->Ax
  //         w->b  contains all of d->b
  //         w->c  contains all of d->c
  idxint i, j, ind, cone, cone_idx;
  double max_val;
  idxint Anz = d->Ap[d->n];

  for(i = 0; i < d->n; ++i) { // cols
    ind = k->f + k->l;
    j = d->Ap[i];

    // seek to the first cone
    while( j < d->Ap[i+1] && d->Ai[j] < ind ) {
      j++;
    }

    cone = 0;
    while (j < d->Ap[i+1]) {
      // while we're in a cone
      max_val = 0;

      // find the maximum in this column
      cone_idx = j;
      while ( cone_idx < Anz &&
              ind <= d->Ai[cone_idx] && 
              d->Ai[cone_idx] < ind + k->q[cone]) 
      {
        max_val = MAX(max_val, fabs(w->Ax[cone_idx]));
        cone_idx++;
      }
      // set all the nonzeros to the maximum value
      // also set all the rows (Ai) to be at ind, that is, this row contains
      // everything we need, the other n-1 rows are just 0.
      // this makes duplicates, but doesn't matter
      while ( j < Anz && 
              ind <= d->Ai[j] && 
              d->Ai[j] < ind + k->q[cone]) 
      {
        w->Ax[j] = max_val;
        w->Ai[j] = ind;
        j++;
      }

      ind += k->q[cone++];
    } // otherwise, move on to next column
  }

  // set b in the same way
  ind = k->f + k->l;
  for(cone = 0; cone < k->qsize; ++cone) {
    max_val = 0;
    for(i = 0; i < k->q[cone]; ++i) {
      max_val = MAX(max_val, fabs(w->b[ind + i]));
    }
    for(i = 0; i < k->q[cone]; ++i) {
      w->b[ind + i] = max_val;
    }
    ind += k->q[cone];
  }
}

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
  w->stilde = w->x + d->n;
  // allocate workspace memory for s and y
  w->s = PDOS_calloc(d->m,sizeof(double));
  w->y = PDOS_calloc(d->m,sizeof(double));

  // allocate workspace memory for normalization matrices
  w->D = PDOS_malloc(d->m*sizeof(double));
  w->E = PDOS_malloc(d->n*sizeof(double));

  if(d->p->NORMALIZE) {
    idxint iters = 0;

    w->Ax = PDOS_calloc(Anz, sizeof(double));
    // create a temporary copy of row indices
    w->Ai = PDOS_calloc(Anz, sizeof(idxint));
    memcpy(w->Ai, d->Ai, Anz*sizeof(idxint));

    w->Ap = d->Ap;
    w->b = PDOS_calloc(d->m, sizeof(double));
    w->c = PDOS_calloc(d->n, sizeof(double));

    double *pi = PDOS_calloc(d->m, sizeof(double));
    for( i=0; i < d->m; ++i ) w->D[i] = 1.0;

    double *delta = PDOS_calloc(d->n, sizeof(double));
    for( i=0; i < d->n; ++i ) w->E[i] = 1.0;

    // since we equilibrate [A b; c^T 0], these correspond to the last row
    // and last column, respectively. we discard them later.
    double lastD = 1.0, lastE = 1.0;
    double lastPi = 0.0, lastDelta = 0.0;

    // set w->Ax = d->Ax
    memcpy(w->Ax, d->Ax, Anz*sizeof(double));

    // set w->c = d->c
    memcpy(w->c, d->c, d->n*sizeof(double));

    // set w->b = d->b
    memcpy(w->b, d->b, d->m*sizeof(double));

    // replace all nonzero entires of w->Ax with max of absolute value of each cone, same with w->b
    // this preserves the cone boundaries
    if (k->qsize > 0) {
      collapseWorkspaceData(w,d,k);
    }

    // for (i = 0; i < Anz; ++i) {
    //   printf("%f ", w->Ax[i]);
    // }
    // printf("\n");

    idxint ind, cone;

    // equilibrate [A b; c^T 0]
    for( iters = 0; iters < EQUILIBRATE_ITERS; ++iters) {

      // compute max across rows
      for(i = 0; i < d->n; ++i) { // cols
        for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
          // IMPORTANT: must use w->Ax and w->Ai since that has collapsed A
          pi[w->Ai[j]] = MAX(pi[w->Ai[j]], fabs(w->Ax[j])*w->E[i]);
        }
        // last row
        lastPi = MAX(lastPi, fabs(w->c[i])*w->E[i]);
      }
      // walk through d->m rows of pi and update with b
      for(i = 0; i < d->m; ++i) { // rows
        pi[i] = MAX(pi[i], fabs(w->b[i])*lastE);
      }

      // now collapse cones together
      // recall that all other rows were just 0
      ind = k->f + k->l;
      for(cone = 0; cone < k->qsize; ++cone) {
        for(i = 0; i < k->q[cone]; ++i) {
          pi[ind+i] = pi[ind];
        }
        ind += k->q[cone];
      }

      for(i = 0; i < d->m; ++i) {
        w->D[i] = sqrt(w->D[i] / pi[i]);
        pi[i] = 0.0;  // set to 0 to compute max
      }
      // handle last element
      lastD = sqrt( lastD / lastPi );
      lastPi = 0.0;

      // now compute max down through columns
      for(i = 0; i < d->n; ++i) { // cols
        for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
          // IMPORTANT: must use w->Ax and w->Ai, since that has collapsed A
          delta[i] = MAX(delta[i], fabs(w->Ax[j])*w->D[w->Ai[j]]);
        }
        // walk through d->n columns of delta and update with c
        delta[i] = MAX(delta[i], fabs(w->c[i])*lastD);
      }
      // last column
      for(i = 0; i < d->m; ++i) { // rows
        lastDelta = MAX(lastDelta, fabs(w->b[i])*w->D[i]);
      }

      for(i = 0; i < d->n; ++i) {
        w->E[i] = sqrt(w->E[i] / delta[i]);
        delta[i] = 0.0;
      }
      // handle last element
      lastE = sqrt( lastE / lastDelta );
      lastDelta = 0.0;
    }

    PDOS_free(pi); PDOS_free(delta);
    // free the temporary copy
    PDOS_free(w->Ai);
    // now set w->Ai properly
    w->Ai = d->Ai;

    // now scale A
    for(i = 0; i < d->n; ++i) { // cols
      for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
        // IMPORTANT: must use d->Ax since w->Ax contains collapsed data
        w->Ax[j] = w->D[d->Ai[j]]*d->Ax[j]*w->E[i];
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

  printf("||b||_2: %f ||c||_2: %f\n", calcNorm(w->b,w->m), calcNorm(w->c,w->n));
  w->lambda = calcNorm(w->b,w->m) / calcNorm(w->c,w->n) ;

  // set ratio of "x" space penalty (1e-6) to "s,y" space penatly (1)
  for( i=0; i < d->n; ++i ) {
    w->E[i] *= SQRT_RATIO;
    w->c[i] *= SQRT_RATIO;
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
