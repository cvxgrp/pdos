#ifndef LINALG_H_GUARD
#define LINALG_H_GUARD

#include <math.h>
#include "pdos.h"

/*
 * All basic linear operations are inlined (and further optimized) by the 
 * compiler. If compiling without optimization, causes code bloat.
 */

// x = b*a
static inline void setAsScaledArray(double *x, const double * a,const double b,idxint len) {
  idxint i;
  for( i=0;i<len;++i ) x[i] = b*a[i];
}

// a*= b
static inline void scaleArray(double * a,const double b,idxint len){
  idxint i;
  for( i=0;i<len;++i) a[i]*=b;
}

// x'*y
static inline double innerProd(const double * x, const double * y, idxint len){
  idxint i;
  double ip = 0.0;
  for ( i=0;i<len;++i){
    ip += x[i]*y[i];
  }
  return ip;
}

// ||v||_2^2
static inline double calcNormSq(const double * v,idxint len){
  idxint i;
  double nmsq = 0.0;
  for ( i=0;i<len;++i){
    nmsq += v[i]*v[i];
  }
  return nmsq;
}

// ||v||_2
static inline double calcNorm(const double * v,idxint len){
  return sqrt(calcNormSq(v, len));
}

// ||v||_inf
static inline double calcNormInf(const double *v, idxint len) {
  double value, max = 0;
  idxint i;
  // compute norm_inf
  for(i = 0; i < len; ++i) {
    value = fabs(v[i]);
    max = (value > max) ? value : max;
  }
  return max;
}

// saxpy a += sc*b
static inline void addScaledArray(double * a, const double * b, idxint n, const double sc){
  idxint i;
  for (i=0;i<n;++i){
    a[i] += sc*b[i];
  }
}

// y += alpha*A*x
static inline void accumByScaledA(const Data *d, const double *x, const double sc, double *y){
  // assumes memory storage exists for y
  
  /* y += A*x */
  idxint p, j, n, *Ap, *Ai ;
  double *Ax ;
  n = d->n ; Ap = d->Ap ; Ai = d->Ai ; Ax = d->Ax ;

  idxint c1, c2;
  
  for (j = 0 ; j < n ; j++)
  {
    c1 = Ap[j]; c2 = Ap[j+1];
    for (p = c1 ; p < c2 ; p++)        
    {   
      y[Ai[p]] += sc * Ax[p] * x[ j ] ;
    }
  }
}

// y += alpha*A'*x
static inline void accumByScaledATrans(const Data *d, const double *x, const double sc, double *y){
  // assumes memory storage exists for y
  
  /* y += A'*x */
  idxint p, j, n, *Ap, *Ai ;
  double *Ax ;
  n = d->n ; Ap = d->Ap ; Ai = d->Ai ; Ax = d->Ax ;

  idxint c1, c2;
  
  for (j = 0 ; j < n ; j++)
  {
    c1 = Ap[j]; c2 = Ap[j+1];
    for (p = c1 ; p < c2 ; p++)        
    {   
      y[j] += sc * Ax[p] * x[ Ai[p] ] ;
    }
  }
}

// y = A*x
static inline void multByA(const Data *d, const double *x, double *y){
  // assumes memory storage exists for y
  
  /* y = A*x */
  idxint p, j, n, *Ap, *Ai ;
  double *Ax ;
  n = d->n ; Ap = d->Ap ; Ai = d->Ai ; Ax = d->Ax ;

  idxint c1, c2;
  
  memset(y,0,d->m*sizeof(double));
  
  for (j = 0 ; j < n ; j++)
  {
    c1 = Ap[j]; c2 = Ap[j+1];
    for (p = c1 ; p < c2 ; p++)        
    {   
      y[Ai[p]] += Ax[p] * x[ j ] ;
    }
  }
}

// y += A*x
static inline void accumByA(const Data *d, const double *x, double *y) {
  accumByScaledA(d,x,1.0,y);
}

// y += A'*x
static inline void accumByATrans(const Data *d, const double *x, double *y) {
  accumByScaledATrans(d,x,1.0,y);
}

// y -= A*x
static inline void decumByA(const Data *d, const double *x, double *y) {
  accumByScaledA(d,x,-1.0,y);
}

// y -= A'*x
static inline void decumByATrans(const Data *d, const double *x, double *y) {
  accumByScaledATrans(d,x,-1.0,y);
}

// norm(A*x + s - b, 'inf')/normA
static inline double calcPriResid(const Data *d, Work *w) {
  // idxint i = 0;
  // for(i = 0; i < d->m; ++i) {
  //   // using stilde as temp vector
  //   w->stilde[i] = w->s[i] - d->b[i];
  // }
  // accumByA(d, w->x, w->stilde);
  
  // equiv to -(A*x + s - b)
  addScaledArray(w->stilde, w->s, d->m, -1); 

  return calcNormInf(w->stilde, d->m); // TODO: normalize by "normA"
}

// norm(A*y + c, 'inf')/normB
static inline double calcDualResid(const Data *d, Work *w) {
  //
  // WARNING: this function must be called *after* calcPriResid
  //
  // assumes stilde allocates max(d->m,d->n) memory
  memcpy(w->stilde, d->c, (d->n)*sizeof(double));
  accumByATrans(d, w->y, w->stilde);
  return calcNormInf(w->stilde, d->n); // TODO: normalize by "normB"
}

// c'*x
static inline double calcPriObj(const Data *d, Work *w) {
  return innerProd(d->c, w->x, d->n);
}

// -b'*y
static inline double calcDualObj(const Data *d, Work *w) {
  return -innerProd(d->b, w->y, d->m);
}


// // x = b*a
// void setAsScaledArray(double *x, const double * a,const double b,idxint len);
// 
// // a*= b
// void scaleArray(double * a,const double b,idxint len);
// 
// // x'*y
// double innerProd(const double * x, const double * y, idxint len);
// 
// // ||v||_2^2
// double calcNormSq(const double * v,idxint len);
// 
// // ||v||_2
// double calcNorm(const double * v,idxint len);
// 
// // ||v||_inf
// double calcNormInf(const double *v, idxint len);
// // saxpy
// void addScaledArray(double * a, const double * b, idxint n, const double sc);
// // y += A*x
// void accumByA(const Data *d, const double *x, double *y);
// // y += A'*x
// void accumByATrans(const Data *d, const double *x, double *y);
// // y -= A*x
// void decumByA(const Data *d, const double *x, double *y);
// // y -= A'*x
// void decumByATrans(const Data *d, const double *x, double *y);
// 
// // norm(A*x + s - b, 'inf')/normA
// double calcPriResid(const Data *d, Work *w);
// // norm(-A*y - c, 'inf')/normB
// double calcDualResid(const Data *d, Work *w);
// // c'*x + b'*y
// double calcSurrogateGap(const Data *d, Work *w);
// 
// double calcCertPriResid(const Data *d, Work *w);
// double calcCertDualResid(const Data *d, Work *w);
// double calcCertPriObj(const Data *d, Work *w);
// double calcCertDualObj(const Data *d, Work *w);
// 

#endif
