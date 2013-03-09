#ifndef coneOS_H_GUARD                                                              
#define coneOS_H_GUARD

/* struct that containing standard problem data */
typedef struct PROBLEM_DATA {
  int n, m; /* problem dimensions */
  /* problem data, A, b, c: */
  double * Ax;
  int * Ai, * Ap;
  double * b, * c;
  int MAX_ITERS, CG_MAX_ITS;
  double EPS_ABS, ALPH, CG_TOL;
  bool VERBOSE;
} Data;

typedef struct SOL_VARS {
  double * x, * z;
  char * status;
} Sol;

typedef struct PRIVATE_DATA Priv;

typedef struct WORK {
  int l;
  // primal variables: z_half and z, dual variable: u
  double *z_half, *z, *u;
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
  
  Priv * p;
} Work;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>    
#include <sys/time.h>
#include <math.h>
#include "cones.h"
#include "util.h"
#include "linAlg.h"

Work * initWork(Data * d);
double calcPriResid(Data * d, Work * w);
void getSolution(Data* d,Work * w,Sol* sol);
void relax(Data * d, Work * w);
Sol * coneOS(Data * d, Cone * k);
void updateDualVars(Work * w);
void projectCones(Data * d,Work * w,Cone * k);
void projectLinSys(Data * d, Work * w);
void setz(Data * d, Work * w, Sol * sol);
void setx(Data * d, Work * w, Sol * sol);
void printSummary(Data * d,Work * w,int i, double err, double EPS_PRI);
void printSol(Data * d, Sol * sol);
void freeWork(Work * w);
void freePriv(Work * w);
#endif
