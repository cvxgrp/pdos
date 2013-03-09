#ifndef PRIV_H_GUARD                                                              
#define PRIV_H_GUARD

#include "cs.h"
#include "coneOS.h"

struct PRIVATE_DATA{
	cs * Q;
  // XXX: add state for lambda
};

Work * initWork(Data* d);
void formQ(Data * d, Work * w);
void freePriv(Work * w);
void projectLinSys(Data * d,Work * w);
void cgCustom(cs *Q,double * b,double * x,int max_its,double tol);
void multByQ(const cs *Q, const double *x, double *y);

#endif
