#ifndef PRIV_H_GUARD                                                              
#define PRIV_H_GUARD

#include "cs.h"
#include "amd.h"
#include "ldl.h"
#include "coneOS.h"

struct PRIVATE_DATA {
	cs * L; /* KKT, and factorization matrix L resp. */
	double * D; /* diagonal matrix of factorization */
	int * P; /* permutation of KKT matrix for factorization */

  // precompue A*c and A'*b
  double *Ac, *Atb;

  // stuff for schur complements
  double *alpha1, *alpha2;
  double s;
};

// XXX: should be named LDL
// also, these routines don't need to be "public"
//void choleskyInit(cs * A, int P[], double **info);
//void choleskyFactor(cs * A, int P[], int Pinv[], cs ** L, double **D);
//void choleskySolve(double *x, double b[], cs * L, double D[], int P[]);
Work * initWork(Data * d);
void freePriv(Work * w);
//cs * formKKT(Data * d, Work * w);
//void factorize(Data * d,Work * w);
#endif
