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

inline Work *commonWorkInit(Data *d) {
  int n_plus_m = d->n + d->m;
  Work * w = PDOS_malloc(sizeof(Work));
  w->l = 2*n_plus_m;
  w->si = d->n;
  w->ri = n_plus_m;
  w->yi = n_plus_m + d->n;
  w->z_half = PDOS_malloc(sizeof(double)*w->l);
  w->z = PDOS_calloc(w->l,sizeof(double));
  w->u = PDOS_calloc(w->l,sizeof(double));  
  w->ztmp = PDOS_calloc(w->l,sizeof(double));

  if(d->NORMALIZE) {
    int i;
    int Anz = d->Ap[d->n];
    // scale A,b,c
    double ds, ps, normA = 0.0;
    // frobenius norm
    for(i = 0; i < Anz; ++i) {
      normA = (normA > fabs(d->Ax[i])) ? normA : fabs(d->Ax[i]);
      //normA += (d->Ax[i]*d->Ax[i]);//((double)d->m*d->n);
    }
    normA = sqrt(normA);
    ds = pow((double)1.0/normA, (double)(d->n)/((double)(d->m + d->n)));
    ps = pow((double)1.0/normA, (double)(d->m)/((double)(d->m + d->n)));

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
