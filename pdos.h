#ifndef PDOS_H_GUARD
#define PDOS_H_GUARD

// include different things depending on which version to compile instead of
// using a different external library
// #include direct/private.h
// #include indirect/private.h

// use the same type for index representation as in SuiteSparse
// coincidentally, Python uses the same type (long)
#include "globals.h"

/* struct containing algorithm parameters */
typedef struct PROBLEM_PARAMS {
  idxint MAX_ITERS, CG_MAX_ITS;
  double EPS_ABS, ALPHA, CG_TOL;
  idxint VERBOSE, NORMALIZE;  // boolean
} Params;

/* struct containing standard problem data */
typedef struct PROBLEM_DATA {
  idxint n, m; /* problem dimensions */
  /* problem data, A, b, c: */
  double * Ax;
  idxint * Ai, * Ap;
  double * b, * c;

  Params * p;
} Data;

typedef struct SOL_VARS {
  idxint n, m; /* solution dimensions */
  double *x, *s, *y;
  char status[16];
} Sol;

typedef struct PRIVATE_DATA Priv;

// contains all the data for the solver
typedef struct WORK {
  idxint n, m; /* problem dimensions */

  // problem data
  // if normalized, this is D*A*E, otherwise just points to data's Ax
  double * Ax;
  idxint * Ai, * Ap;  // these just point to data's Ai and Ap
  // stores the transpose of the sparse matrix A
  double * Atx;
  idxint * Ati, * Atp;

  // b, c
  double * b, * c;    // if normalized, these are D*b, E*c


  // x in R^n
  // s, y, stilde in R^m
  // x points to a vector that is n+m long, stilde points to the bottom 'm'
  // of that vector
  double *x, *s, *y, *stilde;

  double *D, *E;  // A = D*A*E

  double lambda;

  Params *params;
  Priv * p;
} Work;


#include <string.h>
#include <sys/time.h>
#include <math.h>
#include "cones.h"
#include "util.h"
#include "linAlg.h"
#include "common.h"

// these are actually library "api"'s
Sol * pdos(const Data * d, const Cone * k);
void freeData(Data **d, Cone **k);
void freeSol(Sol **sol);

// these are pulled in from private.o
Work * initWork(const Data * d, const Cone * k);
void projectLinSys(Work * w);
void freePriv(Work *w);


#endif
