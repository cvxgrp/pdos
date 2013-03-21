#include "linAlg.h"

// x = b*a
extern void setAsScaledArray(double *x, const double * a,const double b,int len);

// a*= b
extern void scaleArray(double * a,const double b,int len);

// x'*y
extern double innerProd(const double * x, const double * y, int len);

// ||v||_2^2
extern double calcNormSq(const double * v,int len);

// ||v||_2
extern double calcNorm(const double * v,int len);

// ||v||_inf
extern double calcNormInf(const double *v, int len);

// saxpy
extern void addScaledArray(double * a, const double * b, int n, const double sc);

// y += alpha*A*x
extern void accumByScaledA(const Data *d, const double *x, const double sc, double *y);

// y += alpha*A'*x
extern void accumByScaledATrans(const Data *d, const double *x, const double sc, double *y);

extern void accumByA(const Data *d, const double *x, double *y);

extern void accumByATrans(const Data *d, const double *x, double *y);

extern void decumByA(const Data *d, const double *x, double *y);

extern void decumByATrans(const Data *d, const double *x, double *y);

// norm(A*x + s - b, 'inf')/normA
extern double calcPriResid(const Data *d, Work *w);

// norm(A*y + c, 'inf')/normB
extern double calcDualResid(const Data *d, Work *w);
// c'*x + b'*y
extern double calcSurrogateGap(const Data *d, Work *w);

extern double calcCertPriResid(const Data *d, Work *w);

// norm(A'*w, 'inf')/normB
extern double calcCertDualResid(const Data *d, Work *w);

// c'*u
extern double calcCertPriObj(const Data *d, Work *w);
// -b'*w
extern double calcCertDualObj(const Data *d, Work *w);
