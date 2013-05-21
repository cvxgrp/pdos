#ifndef PRIV_H_GUARD                                                              
#define PRIV_H_GUARD

#include "cs.h"
#include "amd.h"
#include "ldl.h"
#include "pdos.h"

struct PRIVATE_DATA {
	cs * L; /* KKT, and factorization matrix L resp. */
	double * D; /* diagonal matrix of factorization */
	idxint * P; /* permutation of KKT matrix for factorization */
};

// XXX: should be named LDL
// also, these routines don't need to be "public"
//void choleskyInit(cs * A, idxint P[], double **info);
//void choleskyFactor(cs * A, idxint P[], idxint Pinv[], cs ** L, double **D);
//void choleskySolve(double *x, double b[], cs * L, double D[], idxint P[]);
Work * initWork(const Data * d);
void freePriv(Work * w);
void projectLinSys(Work * w);
//cs * formKKT(Data * d, Work * w);
//void factorize(Data * d,Work * w);


#endif
