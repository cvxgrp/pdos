#ifndef LINALG_H_GUARD
#define LINALG_H_GUARD

#include <math.h>
#include "coneOS.h"

// x = b*a
void setAsScaledArray(double *x, const double * a,const double b,int len);
// a *= b
void scaleArray(double * a,const double b,int len);
// x'*y
double innerProd(const double * x, const double * y, int len);
double calcNorm(const double * v,int len);
double calcNormSq(const double *v, int len);
double calcPriResid(Data* d,Work * w);
// a += sc*b
void addScaledArray(double * a, const double * b, int n, const double sc);
// y += A*x
void accumByA(const Data *d, const double *x, double *y);
// y += A'*x
void accumByATrans(const Data *d, const double *x, double *y);
// y -= A*x
void decumByA(const Data *d, const double *x, double *y);
// y -= A'*x
void decumByATrans(const Data *d, const double *x, double *y);

#endif
