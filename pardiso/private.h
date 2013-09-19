#ifndef PRIV_H_GUARD                                                              
#define PRIV_H_GUARD

//#include "cs.h"
#include "pdos.h"

/* private data for Pardiso solver */
struct PRIVATE_DATA{
  /* Internal solver memory pointer pt,                  */
  /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
  /* or void *pt[64] should be OK on both architectures  */ 
  void    *pt[64]; 

  /* Pardiso control parameters. */
  int      iparm[64];
  double   dparm[64];
  int      maxfct, mnum, phase, error, msglvl, solver;

  /* Matrix type */
  int    mtype;    

  /* Number of equations */
  int    n;

  /* Data for [-I A'; A I] system in *row* compressed form.*/
  /* Should only contain upper triangular component. */
  /* Also needs to be 1-based indexing. */
  int    *Ui;  /* row pointers */
  int    *Uj;  /* column values */
  double *U;
  int     Unz; /* number of nonzeros */
  double *tmp; /* temporary workspace for pardsio? */
};

Work * initWork(const Data* d, const Cone *k);
void freePriv(Work * w);
void projectLinSys(Work * w);

#endif
