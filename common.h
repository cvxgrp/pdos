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
#define MAX(a,b) (a) > (b) ? (a) : (b)
#endif

static inline Work *commonWorkInit(const Data *d) {
  Work * w = PDOS_malloc(sizeof(Work));
  // ensure that x, s are contiguous in memory
  w->x = PDOS_calloc(d->n + d->m,sizeof(double));
  w->s = w->x + d->n;
  w->stilde = PDOS_calloc(MAX(d->m,d->n),sizeof(double));
  w->y = PDOS_calloc(d->m,sizeof(double));  
  //w->xtmp = PDOS_calloc(d->n,sizeof(double));

  if(d->NORMALIZE) {
    idxint i,j,k = 0;
    idxint Anz = d->Ap[d->n];
    double *rowsum = PDOS_calloc(d->m, sizeof(double));
    double *colsum = PDOS_calloc(d->n, sizeof(double));
    // scale A,b,c
    double ds, ps, normA = 0.0,normB = 0.0;
    
    for(i = 0; i < d->n; ++i) { // cols
      for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
        rowsum[d->Ai[k]] += fabs(d->Ax[k]);
        colsum[i] += fabs(d->Ax[k]);            
        k++;
      }
    }
    
    // normA is max column sum
    for(i = 0; i < d->n; ++i) {
      normA = (normA > colsum[i]) ? normA : colsum[i];
    }
    
    // normB is max row sum
    for(i = 0; i < d->m; ++i) {
      normB = (normB > rowsum[i]) ? normB : rowsum[i];
    }
    
    PDOS_free(rowsum); PDOS_free(colsum);
    
    ds = pow((double)d->n/normA, (double)(d->n)/((double)(d->m + d->n)));
    ps = pow((double)d->m/normB, (double)(d->m)/((double)(d->m + d->n)));

    for(i = 0; i < Anz; ++i) {
      d->Ax[i] *= ds*ps;
    }
    for(i = 0; i < d->n; ++i) {
      d->c[i] *= ps;
    }
    for(i = 0; i < d->m; ++i) {
      d->b[i] *= ds;
    }
    w->dual_scale = ds;
    w->primal_scale = ps;
  
  } else {
    w->dual_scale = 1.0;
    w->primal_scale = 1.0;
  }
  return w;
}

#endif
