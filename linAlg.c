#include "linAlg.h"

void scaleArray(double * a,double b,int len){
	int i;
	for( i=0;i<len;i++) a[i]*=b;
}

double innerProd(double * x, double * y, int len){
	int i;
	double ip = 0.0;
	for ( i=0;i<len;i++){
		ip += x[i]*y[i];
	}
	return ip;
}

double calcNorm(double * v,int len){
	int i;
	double nm = 0.0;
	for ( i=0;i<len;i++){
		nm+=v[i]*v[i];
	}
	return sqrt(nm);
}

double calcPriResid(Data* d,Work * w){
	double err = 0.0, tmp;
	int i;
	for ( i=0;i<w->l;i++){
		tmp = (w->uv[i] - w->uv_t[i]);
		err += tmp * tmp;
	}
	return sqrt(err);
}

void addScaledArray(double * a, const double * b, int n, double sc){
	int i;
	for (i=0;i<n;i++){
		a[i] += sc*b[i];
	}
}
