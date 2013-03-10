#ifndef LINALG_H_GUARD
#define LINALG_H_GUARD

#include <math.h>
#include "pdos.h"

/*
 * All basic linear operations are inlined (and further optimized) by the 
 * compiler. If compiling without optimization, causes code bloat.
 */

// x = b*a
static inline void setAsScaledArray(double *x, const double * a,const double b,int len) {
  int i;
	for( i=0;i<len;++i ) x[i] = b*a[i];
}

// a*= b
static inline void scaleArray(double * a,const double b,int len){
	int i;
	for( i=0;i<len;++i) a[i]*=b;
}

// x'*y
static inline double innerProd(const double * x, const double * y, int len){
	int i;
	double ip = 0.0;
	for ( i=0;i<len;++i){
		ip += x[i]*y[i];
	}
	return ip;
}

// ||v||_2^2
static inline double calcNormSq(const double * v,int len){
	int i;
	double nmsq = 0.0;
	for ( i=0;i<len;++i){
		nmsq += v[i]*v[i];
	}
	return nmsq;
}

// ||v||_2
static inline double calcNorm(const double * v,int len){
	return sqrt(calcNormSq(v, len));
}

// ||v||_inf
static inline double calcNormInf(const double *v, int len) {
  double value, max = 0;
  int i;
  // compute norm_inf
  for(i = 0; i < len; ++i) {
    value = abs(v[i]);
    max = (value > max) ? value : max;
  }
  return max;
}

// saxpy
static inline void addScaledArray(double * a, const double * b, int n, const double sc){
	int i;
	for (i=0;i<n;++i){
		a[i] += sc*b[i];
	}
}

// y += alpha*A*x
static inline void accumByScaledA(const Data *d, const double *x, const int sc, double *y){
  // assumes memory storage exists for y
  
	/* y += A*x */
	int p, j, n, m, *Ap, *Ai ;
	double *Ax ;
	m = d->m; n = d->n ; Ap = d->Ap ; Ai = d->Ai ; Ax = d->Ax ;

	int c1, c2;
  
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
static inline void accumByScaledATrans(const Data *d, const double *x, const int sc, double *y){
  // assumes memory storage exists for y
  
	/* y += A'*x */
	int p, j, n, *Ap, *Ai ;
	double *Ax ;
	n = d->n ; Ap = d->Ap ; Ai = d->Ai ; Ax = d->Ax ;

	int c1, c2;
  
	for (j = 0 ; j < n ; j++)
	{
		c1 = Ap[j]; c2 = Ap[j+1];
		for (p = c1 ; p < c2 ; p++)        
		{   
			y[j] += sc*Ax[p] * x[ Ai[p] ] ;
		}
	}
}

static inline void accumByA(const Data *d, const double *x, double *y) {
  accumByScaledA(d,x,1,y);
}

static inline void accumByATrans(const Data *d, const double *x, double *y) {
  accumByScaledATrans(d,x,1,y);
}

static inline void decumByA(const Data *d, const double *x, double *y) {
  accumByScaledA(d,x,-1,y);
}

static inline void decumByATrans(const Data *d, const double *x, double *y) {
  accumByScaledATrans(d,x,-1,y);
}

// norm(A*x + s - b, 'inf')/normA
static inline double calcPriResid(const Data *d, Work *w) {
  int i = 0;
  for(i = w->si; i < w->ri; ++i) {
    w->ztmp[i] = w->z[i] - d->b[i - w->si];
  }
  accumByA(d, w->z, w->ztmp + (w->si));

  return calcNormInf(w->ztmp + (w->si), d->m); // TODO: normalize by "normA"
}

// norm(A*y + c, 'inf')/normB
static inline double calcDualResid(const Data *d, Work *w) {
  memcpy(w->ztmp, d->c, (d->n)*sizeof(double));
  accumByATrans(d, w->z + (w->yi), w->ztmp);
  return calcNormInf(w->ztmp, d->n); // TODO: normalize by "normB"
}

// c'*x + b'*y
static inline double calcSurrogateGap(const Data *d, Work *w) {
  return innerProd(d->c, w->z, d->n) + innerProd(d->b, w->z + (w->yi), d->m);
}


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



#endif
