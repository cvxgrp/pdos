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
};

void choleskyInit(int n, cs * A, int P[], double **info);
void choleskyFactor(int n, cs * A, int P[], int Pinv[], cs ** L, double **D);
void choleskySolve(int n, double *x, double b[], cs * L, double D[], int P[]);
Work * initWork(Data * d);
void freePriv(Work * w);
cs * formKKT(Data * d, Work * w);
void factorize(Data * d,Work * w);
#endif
