#ifndef COMMON_H_GUARD                                                              
#define COMMON_H_GUARD

// redefine printfs and memory allocators as needed
#ifdef MATLAB_MEX_FILE
  #include "mex.h"
  #define PDOS_printf   mexPrintf
  #define PDOS_free     mxFree
  #define PDOS_malloc   mxMalloc
  #define PDOS_calloc   mxCalloc
#elif defined PYTHON
  #include <Python.h>
  #include <stdlib.h>
  #define PDOS_printf   PySys_WriteStdout
  #define PDOS_free     free
  #define PDOS_malloc   malloc
  #define PDOS_calloc   calloc
#else
  #include <stdio.h>
  #include <stdlib.h>
  #define PDOS_printf   printf
  #define PDOS_free     free
  #define PDOS_malloc   malloc
  #define PDOS_calloc   calloc
#endif

#include "pdos.h"
#include <math.h>

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#define EQUILIBRATE_ITERS 3

static inline Work *commonWorkInit(const Data *d, const Cone *k) {
  Work * w = PDOS_malloc(sizeof(Work));
  // copy dimensions
  w->m = d->m; w->n = d->n;
  // copy parameters (pointer)
  w->params = d->p;
  // ensure that x, stilde are contiguous in memory
  w->x = PDOS_calloc(d->n + MAX(d->m,d->n),sizeof(double));
  w->stilde = w->x + d->n;
  w->s = PDOS_calloc(d->m,sizeof(double));
  w->y = PDOS_calloc(d->m,sizeof(double));  

  if(d->p->NORMALIZE) {
    idxint i,j = 0;
    idxint iters = 0;
    idxint Anz = d->Ap[d->n];
    
    w->Ax = PDOS_calloc(Anz, sizeof(double));
    w->Ai = d->Ai;
    w->Ap = d->Ap;
    w->b = PDOS_calloc(d->m, sizeof(double));
    w->c = PDOS_calloc(d->n, sizeof(double));
    
    w->D = PDOS_malloc(d->m*sizeof(double));
    double *pi = PDOS_calloc(d->m, sizeof(double));
    for( i=0; i < d->m; ++i ) w->D[i] = 1.0;

    w->E = PDOS_malloc(d->n*sizeof(double));
    double *delta = PDOS_calloc(d->n, sizeof(double));
    for( i=0; i < d->n; ++i ) w->E[i] = 1.0;
    
    // set w->Ax = d->Ax
    for(i = 0; i < Anz; ++i) w->Ax[i] = d->Ax[i];
    
    
    for( iters = 0; iters < EQUILIBRATE_ITERS; ++iters) {
      // // As = As*diag(delta)
      // for(i = 0; i < d->n; ++i) { // cols
      //   for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
      //     d->Ax[j] *= delta[i];
      //   }
      //   delta[i] = 0.0; // in preparation for computing max
      // }
      // 
      
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
    
    // now compute max down through columns
    for(i = 0; i < d->n; ++i) { // cols
      for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
        w->Ax[j] *= w->D[d->Ai[j]]*w->E[i];            
      }
    }

    for(i = 0; i < d->n; ++i) {
      w->c[i] = d->c[i]*w->E[i];
    }
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
    w->Ax = d->Ax;
    w->Ai = d->Ai;
    w->Ap = d->Ap;
    w->b = d->b;
    w->c = d->c;
    
    idxint i;
    
    w->D = PDOS_malloc(d->m*sizeof(double));
    for( i=0; i < d->m; ++i ) w->D[i] = 1.0;

    w->E = PDOS_malloc(d->n*sizeof(double));
    for( i=0; i < d->n; ++i ) w->E[i] = 1.0;
  }
  
  w->lambda = sqrt( w->n*(1.0 + calcNormSq(w->b,w->m)) / (w->m*(1.0 + calcNormSq(w->c,w->n))) );
    
  return w;
}

#endif
