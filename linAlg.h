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
static inline void accumByScaledA(const Work *w, const double *x, const double sc, double *y){
  // assumes memory storage exists for y
  
  /* y += A*x */
  idxint p, j, n, *Ap, *Ai ;
  double *Ax ;
  n = w->n ; Ap = w->Ap ; Ai = w->Ai ; Ax = w->Ax ;

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
static inline void accumByScaledATrans(const Work *w, const double *x, const double sc, double *y){
  // assumes memory storage exists for y
  
  /* y += A'*x */
  idxint p, j, n, *Ap, *Ai ;
  double *Ax ;
  n = w->n ; Ap = w->Ap ; Ai = w->Ai ; Ax = w->Ax ;

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
static inline void multByA(const Work *w, const double *x, double *y){
  // assumes memory storage exists for y
  
  /* y = A*x */
  idxint p, j, n, m, *Ap, *Ai ;
  double *Ax ;
  n = w->n ; m = w->m; Ap = w->Ap ; Ai = w->Ai ; Ax = w->Ax ;

  idxint c1, c2;
  
  memset(y,0,m*sizeof(double));
  
  for (j = 0 ; j < n ; j++)
  {
    c1 = Ap[j]; c2 = Ap[j+1];
    for (p = c1 ; p < c2 ; p++)        
    {   
      y[Ai[p]] += Ax[p] * x[ j ] ;
    }
  }
}

// y = A'*x
static inline void multByATrans(const Work *w, const double *x, double *y){
  // assumes memory storage exists for y
  
  /* y = A'*x */
  idxint p, j, n, *Ap, *Ai ;
  double *Ax ;
  n = w->n ; Ap = w->Ap ; Ai = w->Ai ; Ax = w->Ax ;

  idxint c1, c2;
  
  for (j = 0 ; j < n ; j++)
  {
    c1 = Ap[j]; c2 = Ap[j+1];
    y[j] = 0.0;
    for (p = c1 ; p < c2 ; p++)        
    {   
      y[j] += Ax[p] * x[ Ai[p] ] ;
    }
  }
}

// y += A*x
static inline void accumByA(const Work *w, const double *x, double *y) {
  accumByScaledA(w,x,1.0,y);
}

// y += A'*x
static inline void accumByATrans(const Work *w, const double *x, double *y) {
  accumByScaledATrans(w,x,1.0,y);
}

// y -= A*x
static inline void decumByA(const Work *w, const double *x, double *y) {
  accumByScaledA(w,x,-1.0,y);
}

// y -= A'*x
static inline void decumByATrans(const Work *w, const double *x, double *y) {
  accumByScaledATrans(w,x,-1.0,y);
}

// norm(A*x + s - b, 'inf')/normA
static inline double calcPriResid(Work *w) {
  idxint i = 0;
  for(i = 0; i < w->m; ++i) {
    // using stilde as temp vector
    w->stilde[i] = w->s[i] - w->b[i];
  }
  accumByA(w, w->x, w->stilde);
  
  // normalize by D
  for(i = 0; i < w->m; ++i) {
    w->stilde[i] /= w->D[i];
  }
  
  // equiv to -(A*x + s - b) when alpha = 1
  //addScaledArray(w->stilde, w->s, d->m, -1); 

  return calcNormInf(w->stilde, w->m);
}

// norm(A*y + c, 'inf')/normB
static inline double calcDualResid(Work *w) {
  // assumes stilde allocates max(d->m,d->n) memory
  memcpy(w->stilde, w->c, (w->n)*sizeof(double));
  accumByATrans(w, w->y, w->stilde);
  
  // normalize by E
  idxint i = 0;
  for(i = 0; i < w->n; ++i) {
    w->stilde[i] /= w->E[i];
  }
  
  return calcNormInf(w->stilde, w->n); // TODO: normalize by "normB"
}

// c'*x
static inline double calcPriObj(const Work *w) {
  return innerProd(w->c, w->x, w->n);
}

// -b'*y
static inline double calcDualObj(const Work *w) {
  return -innerProd(w->b, w->y, w->m);
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
