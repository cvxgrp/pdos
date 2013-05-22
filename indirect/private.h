#ifndef PRIV_H_GUARD                                                              
#define PRIV_H_GUARD

//#include "cs.h"
#include "pdos.h"

struct PRIVATE_DATA{
	//cs * Q; // i don't actually need Q, just need "forward" and "adjoint"
  // XXX: add state for lambda
  // double *lambda;
  double *p;  // cg iterate
  double *q;  // cg residual
  double *Ax; // temp memory for holding result of A*x
};

Work * initWork(const Data* d, const Cone *k);
//void formQ(Data * d, Work * w);
void freePriv(Work * w);
void projectLinSys(Work * w);
//void cgCustom(cs *Q,double * b,double * x,idxint max_its,double tol);
//void multByQ(const cs *Q, const double *x, double *y);

#endif
