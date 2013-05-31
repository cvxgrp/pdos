#ifndef COMMON_H_GUARD                                                              
#define COMMON_H_GUARD

#include "pdos.h"
#include <math.h>

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#define EQUILIBRATE_ITERS 3

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
    idxint i,j = 0;
    idxint iters = 0;
    idxint Anz = d->Ap[d->n];
    
    w->Ax = PDOS_calloc(Anz, sizeof(double));
    w->Ai = d->Ai;
    w->Ap = d->Ap;
    w->b = PDOS_calloc(d->m, sizeof(double));
    w->c = PDOS_calloc(d->n, sizeof(double));
    
    double *pi = PDOS_calloc(d->m, sizeof(double));
    for( i=0; i < d->m; ++i ) w->D[i] = 1.0;

    double *delta = PDOS_calloc(d->n, sizeof(double));
    for( i=0; i < d->n; ++i ) w->E[i] = 1.0;
    
    // set w->Ax = d->Ax
    for(i = 0; i < Anz; ++i) w->Ax[i] = d->Ax[i];
    
    for( iters = 0; iters < EQUILIBRATE_ITERS; ++iters) {
      
      // compute max across rows
      for(i = 0; i < d->n; ++i) { // cols
        for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
          pi[d->Ai[j]] = MAX(pi[d->Ai[j]], fabs(d->Ax[j])*w->E[i]);            
        }
      }
      
      // now collapse cones together
      idxint ind = k->f + k->l;
      for(i = 0; i < k->qsize; ++i) {
        // find the maximum in this cone
        double cone_max = 0.0;
        for(j = ind; j < ind + k->q[i]; ++j) {
          cone_max = MAX(pi[j], cone_max);
        }
        // set all in this cone to the maximum
        for(j = ind; j < ind + k->q[i]; ++j) {
          pi[j] = cone_max;
        }
        ind += k->q[i];
      }
      
      for(i = 0; i < d->m; ++i) {
        w->D[i] = sqrt(w->D[i] / pi[i]);
        pi[i] = 0.0;  // set to 0 to compute max
      }
      
      // now compute max down through columns
      for(i = 0; i < d->n; ++i) { // cols
        for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
          delta[i] = MAX(delta[i], fabs(d->Ax[j])*w->D[d->Ai[j]]);            
        }
      }
      
      for(i = 0; i < d->n; ++i) {
        w->E[i] = sqrt(w->E[i] / delta[i]);
        delta[i] = 0.0;
      }
    }
    
    
    PDOS_free(pi); PDOS_free(delta);
    
    // now scale A
    for(i = 0; i < d->n; ++i) { // cols
      for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
        w->Ax[j] *= w->D[d->Ai[j]]*w->E[i];            
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
    
    // display first col of A
    // for(i = 0; i < d->n; ++i) { // cols
    //   for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
    //     printf("%f ", d->Ax[j]);         
    //   }
    //   printf("\n");
    // }
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
  
  w->lambda = sqrt( w->n*(1.0 + calcNormSq(w->b,w->m)) / (w->m*(1.0 + calcNormSq(w->c,w->n))) );
    
  return w;
}

#endif
