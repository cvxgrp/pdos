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

Work * initWork(const Data * d, const Cone *k);
void freePriv(Work * w);
void projectLinSys(Work * w);

#endif
