#ifndef LINALG_H_GUARD
#define LINALG_H_GUARD

#include <math.h>
#include "coneOS.h"

void scaleArray(double * a,double b,int len);
double innerProd(double * x, double * y, int len);
double calcNorm(double * v,int len);
double calcPriResid(Data* d,Work * w);
void addScaledArray(double * a, const double * b, int n, double sc);

#endif
