#ifndef PRIV_H_GUARD                                                              
#define PRIV_H_GUARD

#include "pdos.h"

struct PRIVATE_DATA{
  double *p;  /* cg iterate */ 
  double *q;  /* cg residual */
  double *Ax; /* temp memory for holding result of A*x */
};

Work * initWork(const Data* d, const Cone *k);
void freePriv(Work * w);
void projectLinSys(Work * w);

#endif
