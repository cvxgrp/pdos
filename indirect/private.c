#include "private.h"

Work * initWork(Data* d){
  Work * w = malloc(sizeof(Work));
  w->l = 2*(d->n + d->m +1);
  w->uv = calloc(w->l,sizeof(double));
  w->uv_t = calloc(w->l,sizeof(double));
  w->lam = calloc(w->l,sizeof(double));
  w->uv_h = malloc(sizeof(double)*w->l);
  w->uv_t[d->n+d->m]=100; w->uv[d->n+d->m]=100;
  w->uv_t[w->l-1]=100; w->uv[w->l-1]=100;
  w->p = malloc(sizeof(Priv));
  formQ(d,w);
  return w;
}

void formQ(Data * d, Work * w){
	int Qnzmax = 2*d->Ap[d->n]+2*d->n+2*d->m;
	cs * Q = cs_spalloc(w->l/2, w->l/2, Qnzmax, 1, 1);
	/* -A in Q : CCS: */
	int j, k, kk=0;
	for (j = 0; j < d->n; j++) {                 
		for (k = d->Ap[j]; k < d->Ap[j+1]; k++) { 
			Q->i[kk] = d->Ai[k] + d->n;
			Q->p[kk] = j;
			Q->x[kk] = -d->Ax[k];
			kk++;
		}   
	}
	/* A^T in Q : CCS: */
	for (j = 0; j < d->n; j++) {                
		for (k = d->Ap[j]; k < d->Ap[j+1]; k++) {
			Q->p[kk] = d->Ai[k] + d->n;
			Q->i[kk] = j;
			Q->x[kk] = d->Ax[k];
			kk++;
		}   
	}
	/* -c^T in Q, same row */
	for(k = 0; k < d->n; k++ ){
		if (d->c[k]!=0){
			Q->i[kk] = d->n + d->m;
			Q->p[kk] = k;
			Q->x[kk] = -d->c[k];
			kk++;
		}
	}
	/* c in Q, same col */
	for(k = 0; k < d->n; k++ ){
		if (d->c[k]!=0){
			Q->i[kk] = k;
			Q->p[kk] = d->n + d->m;
			Q->x[kk] = d->c[k];
			kk++;
		}
	}
	/* -b^T in Q, same row */
	for(k = 0; k < d->m; k++ ){
		if (d->b[k]!=0){
			Q->i[kk] = d->n + d->m;
			Q->p[kk] = k + d->n;
			Q->x[kk] = -d->b[k];
			kk++;
		}
	}
	/* b in Q, same col */
	for(k = 0; k < d->m; k++ ){
		if (d->b[k]!=0){
			Q->i[kk] = k + d->n;
			Q->p[kk] = d->n + d->m;
			Q->x[kk] = d->b[k];
			kk++;
		}
	}
  Q->nz = kk;
	w->p->Q = cs_compress(Q);
	cs_spfree(Q);
}

void freePriv(Work * w){
  cs_spfree(w->p->Q);free(w->p);
}

void projectLinSys(Data * d, Work * w){
	memcpy(w->uv_h,w->uv_t,w->l*sizeof(double));
	addScaledArray(w->uv_h,w->lam,w->l,-1);    

	double * tmp = malloc(sizeof(double)*w->l/2);
	double * f = w->uv_h, *g = &w->uv_h[w->l/2];
	multByQ(w->p->Q,g,tmp);
	addScaledArray(f,tmp,w->l/2,-1);
	/* warm start cg with uv_t */
	memcpy(w->uv,w->uv_t,sizeof(double)*w->l);
	cgCustom(w->p->Q, f, w->uv, d->CG_MAX_ITS, d->CG_TOL);
	multByQ(w->p->Q,w->uv,&w->uv[w->l/2]);
	free(tmp);
}     

void cgCustom(cs *Q,double * b,double * x,int max_its,double tol){
	int n = Q->n, i;
	double alpha, rsnew=0;
	double * tmp = malloc(sizeof(double)*n);
	double * p = malloc(sizeof(double)*n);
	double * Ap = malloc(sizeof(double)*n);
	double * r = malloc(sizeof(double)*n);

	multByQ(Q,x,tmp);
	multByQ(Q,tmp,r);
	addScaledArray(r,b,n,1);
	addScaledArray(r,x,n,-1);

	memcpy(p,r,n*sizeof(double));
	double rsold=calcNorm(r,n);
	for (i=0; i< max_its; i++){
		multByQ(Q,p,tmp);
		multByQ(Q,tmp,Ap);
		addScaledArray(Ap,p,n,-1);
		alpha=-(rsold*rsold)/innerProd(p,Ap,n);

		addScaledArray(x,p,n,alpha);
		addScaledArray(r,Ap,n,alpha);    

		rsnew=calcNorm(r,n);
		if (rsnew<tol){
			break;
		}
		scaleArray(p,(rsnew*rsnew)/(rsold*rsold),n);
    addScaledArray(p,r,n,1);
		rsold=rsnew;
	}
	printf("terminating cg residual = %4f, took %i itns\n",rsnew,i);
	free(tmp);free(p);free(Ap);free(r);
}

void multByQ(const cs *Q, const double *x, double *y)
{
	/* y  = Q*x */
	/* parallelizes over rows,
	   relies on special property that Q = -Q^T 
	   we use this so we can parallelize across columns
	   since Q is stored in compressed column format */
	int p, j, n, *Qp, *Qi ;
	double *Qx ;
	n = Q->n ; Qp = Q->p ; Qi = Q->i ; Qx = Q->x ;
	/* parallel matrix multiply */
	int c1, c2;
	double yj;
#pragma omp parallel for private(p,c1,c2,yj) 
	for (j = 0 ; j < n ; j++)
	{
		yj = 0.0;
		c1 = Qp[j]; c2 = Qp[j+1];
		for (p = c1 ; p < c2 ; p++)        
		{   
			yj -= Qx[p] * x[ Qi[p] ] ;
		}
		y[j] = yj;
	}
}
