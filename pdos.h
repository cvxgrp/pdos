#ifndef PDOS_H_GUARD                                                              
#define PDOS_H_GUARD

// include different things depending on which version to compile instead of
// using a different external library
// #include direct/private.h
// #include indirect/private.h

/* struct that containing standard problem data */
typedef struct PROBLEM_DATA {
  int n, m; /* problem dimensions */
  /* problem data, A, b, c: */
  double * Ax;
  int * Ai, * Ap;
  double * b, * c;
  int MAX_ITERS, CG_MAX_ITS;
  double EPS_ABS, ALPH, CG_TOL;
  int VERBOSE;  // boolean
} Data;

typedef struct SOL_VARS {
  double * x, * y;
  char * status;
} Sol;

typedef struct PRIVATE_DATA Priv;

typedef struct WORK {
  int l;
  int si, ri, yi; // xi = 0

  // primal variables: z_half and z, dual variable: u
  // ztmp is just a temporary variable for memory
  double *z_half, *z, *ztmp, *u;
  // pointers to memory locations
  // double *x_half, *s_half, *r_half, *y_half;
  // double *x, *s, *r, *y;
  // double *r_bar, *y_bar, *x_bar, *s_bar;
  Priv * p;
} Work;

// for a first cut, we're just going to use 4*l instead of 2*l storage
// z_half needs "l", z needs "l/2", u needs "l/2"
  
// in matlab notation....
// x_half = x = z_half(1:n)
// s_half = z_half(n+1:n+m)
// r_half = z_half(n+m+1:2*n+m)
// y_half = z_half(2*n+m+1:end)
//
// x = x_half       // can use this space for other stuff (for LDL?)
// s = z(n+1:n+m)
// r = 0            // can use this space for other stuff (for LDL?)
// y = z(2*n+m+1:end)
//
// r_bar = 0
// y_bar = u(n+1:n+m)
// x_bar = u(n+m+1:2*n+m) = running sum(r_half)
// s_bar = u(2*n+m+1:end)
  
// //double *uv,*uv_t,*lam,*uv_h; /* variables */
// int n, m; // do we need it? probably not
// 
// // this assumes that the first n components of the dual variable u
// // are zero to begin with
// //
// // primal variable z = (x,s,r,y)  (with r = 0)
// 
// // we don't represent "x_half" since
// //   x_half = x - r_bar = x (because r_bar = 0)
// double *s_half, *y_half;
// 
// // we don't store "r" since r = 0
// double *x, *s, *y;  // these need memory (x is shared with x_half)
// 
// // x_{k+1/2} = (x_k - r_bar) - A'*d1 + c*d3
// 
// // dual variable u = (r_bar, y_bar, x_bar, s_bar) (with r_bar = 0)
// // the initial dual variable u = 0
// double *y_bar, *x_bar, *s_bar;
// // x_bar = x_bar + r_half

#include <stdio.h>
#include <stdlib.h>
#include <string.h>    
#include <sys/time.h>
#include <math.h>
#include "cones.h"
#include "util.h"
#include "linAlg.h"

//Work * initWork(Data * d);
//double calcPriResid(Data * d, Work * w);
//void getSolution(Data* d,Work * w,Sol* sol);
//void relax(Data * d, Work * w);
Sol * pdos(Data * d, Cone * k);
//void updateDualVars(Work * w);
//void projectCones(Data * d,Work * w,Cone * k);
//void projectLinSys(Data * d, Work * w);
//void sety(Data * d, Work * w, Sol * sol);
//void setx(Data * d, Work * w, Sol * sol);
//void printSummary(Data * d,Work * w,int i, double err, double EPS_PRI);
//void printSol(Data * d, Sol * sol);
//void freeWork(Work * w);
//void freePriv(Work * w);

Work * initWork(Data * d);
void projectLinSys(Data * d, Work * w);
void freePriv(Work * w);
// inline declarations
static inline void relax(Data * d, Work * w);
static inline void updateDualVars(Work * w);
static inline void prepZVariable(Work *w);
static inline void projectCones(Data * d,Work * w,Cone * k);
static inline void sety(Data * d, Work * w, Sol * sol);
static inline void setx(Data * d, Work * w, Sol * sol);
static inline void getSolution(Data * d, Work * w, Sol * sol, int solver_state);
static inline void printSummary(Data * d,Work * w,int i, double p_res, double d_res, double eta);
static inline void printHeader();
static inline void printSol(Data * d, Sol * sol);
static inline void freeWork(Work * w);

#endif
