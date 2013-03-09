#include "private.h"

void freePriv(Work * w){
  cs_spfree(w->p->L);free(w->p->P);free(w->p->D);free(w->p);
}

void projectLinSys(Data *d, Work * w){
  memcpy(w->uv_h,w->uv_t,w->l*sizeof(double));
  addScaledArray(w->uv_h,w->lam,w->l,-1);
  choleskySolve(w->l, w->uv, w->uv_h, w->p->L, w->p->D, w->p->P);
  addScaledArray(&(w->uv[w->l/2]),&(w->uv_h[w->l/2]),w->l/2,1);
}

Work * initWork(Data* d){ 
  Work * w = malloc(sizeof(Work));
  w->p = malloc(sizeof(Priv));
  w->l = 2*(d->n + d->m +1);
  w->uv = malloc(sizeof(double)*w->l);
  w->uv_t = calloc(w->l,sizeof(double));
  w->lam = calloc(w->l,sizeof(double));  
  w->uv_h = malloc(sizeof(double)*w->l);
  w->uv_t[d->n+d->m]=100;
  w->uv_t[w->l-1]=100;

  w->p->P = malloc(sizeof(int)*w->l);
  w->p->L = malloc (sizeof (cs));
  w->p->L->m = w->l;
  w->p->L->n = w->l;
  w->p->L->nz = -1; 
  factorize(d,w);
  return w;
}

cs * formKKT(Data * d, Work * w){
	/* ONLY UPPER TRIANGULAR PART IS STUFFED
	 * forms column compressed KKT matrix
	 * assumes column compressed form A matrix
	 */
	int j, k, kk = 0;
	/* I at top left */
	int Knzmax = w->l+2*d->Ap[d->n]+2*d->n+2*d->m;
	cs * K = cs_spalloc(w->l, w->l, Knzmax, 1, 1);
	for (k = 0; k<w->l/2; k++){
		K->i[kk] = k;
		K->p[kk] = k;
		K->x[kk] = 1;
		kk++;
	}
	/* A in Q^T at top right : CCS: */
	for (j = 0; j < d->n; j++) {                 
		for (k = d->Ap[j]; k < d->Ap[j+1]; k++) { 
			K->i[kk] = d->Ai[k] + d->n;
			K->p[kk] = j + w->l/2;
			K->x[kk] = d->Ax[k];
			kk++;
		}   
	}
	/* if A in TRIPLET:
	   for(k = 0; k < d->A->nzmax; k++ ){
	   K->i[kk] = d->A->i[k] + d->n;
	   K->p[kk] = d->A->p[k] + w->l/2;
	   K->x[kk] = d->A->x[k];
	   kk++;
	   }
	 */
	/* -A^T in Q^T at top right : CCS: */
	for (j = 0; j < d->n; j++) {                
		for (k = d->Ap[j]; k < d->Ap[j+1]; k++) {
			K->p[kk] = d->Ai[k] + d->n + w->l/2;
			K->i[kk] = j;
			K->x[kk] = -d->Ax[k];
			kk++;
		}   
	}
	/* if A in TRIPLET:
	   for(k = 0; k < d->A->nzmax; k++ ){
	   K->i[kk] = d->A->p[k];
	   K->p[kk] = d->A->i[k] + d->n + w->l/2;
	   K->x[kk] = -d->A->x[k];
	   kk++;
	   }
	 */
	/* c^T in Q, same row */
	for(k = 0; k < d->n; k++ ){
		if (d->c[k]!=0){
			K->i[kk] = d->n + d->m;
			K->p[kk] = k + w->l/2;
			K->x[kk] = d->c[k];
			kk++;
		}
	}
	/* -c in Q, same col */
	for(k = 0; k < d->n; k++ ){
		if (d->c[k]!=0){
			K->i[kk] = k;
			K->p[kk] = w->l/2 + d->n + d->m;
			K->x[kk] = -d->c[k];
			kk++;
		}
	}
	/* b^T in Q, same row */
	for(k = 0; k < d->m; k++ ){
		if (d->b[k]!=0){
			K->i[kk] = d->n + d->m;
			K->p[kk] = k + d->n + w->l/2;
			K->x[kk] = d->b[k];
			kk++;
		}
	}
	/* -b in Q, same col */
	for(k = 0; k < d->m; k++ ){
		if (d->b[k]!=0){
			K->i[kk] = k + d->n;
			K->p[kk] = w->l/2 + d->n + d->m;
			K->x[kk] = -d->b[k];
			kk++;
		}
	}
	/* -I at bottom right */
	for (k = 0; k<w->l/2; k++){
		K->i[kk] = k + w->l/2;
		K->p[kk] = k + w->l/2;
		K->x[kk] = -1;
		kk++;
	}
	K->nz = kk;
	cs * K_cs = cs_compress(K);
	cs_spfree(K);
  return(K_cs);
}


void factorize(Data * d,Work * w){
  tic();
  cs * K = formKKT(d,w);
  printf("KKT matrix factorization info:\n");
  double *info;
  choleskyInit(w->l, K, w->p->P, &info);
  amd_info(info);
  int * Pinv = cs_pinv(w->p->P, w->l);
  cs * C = cs_symperm(K, Pinv, 1); 
  choleskyFactor(w->l, C, NULL, NULL, &w->p->L, &w->p->D);
  printf("KKT matrix factorization took %4.2f s\n",tocq());
  cs_spfree(C);cs_spfree(K);free(Pinv);free(info);
}


void choleskyInit(int n, cs * A, int P[], double **info) {
	*info  = (double *) malloc(AMD_INFO * sizeof(double));
	amd_order(n, A->p, A->i, P, (double *) NULL, *info);
}

void choleskyFactor(int n, cs * A, int P[], int Pinv[], cs **L , double **D) 
{
	(*L)->p = (int *) malloc((n+1) * sizeof(int));
	int Parent[n], Lnz[n], Flag[n], Pattern[n];
	double Y[n];

	ldl_symbolic(n, A->p, A->i, (*L)->p, Parent, Lnz, Flag, P, Pinv);

	(*L)->nzmax = *((*L)->p + n);
  (*L)->x = (double *) malloc((*L)->nzmax * sizeof(double));
	(*L)->i =    (int *) malloc((*L)->nzmax * sizeof(int));
  *D  = (double *) malloc(n    * sizeof(double));

	ldl_numeric(n, A->p, A->i, A->x, (*L)->p, Parent, Lnz, (*L)->i, (*L)->x, *D, Y, Pattern, Flag, P, Pinv);
}

void choleskySolve(int n, double *x, double b[], cs * L, double D[], int P[])
{
	double bp[n];
  int i;
	if (P == NULL) {
		for ( i = 0; i < n; i++) { x[i] = b[i]; }
		ldl_lsolve(n, x, L->p, L->i, L->x);
		ldl_dsolve(n, x, D);
		ldl_ltsolve(n, x, L->p, L->i, L->x);
	} else {
		ldl_perm(n, bp, b, P);
		ldl_lsolve(n, bp, L->p, L->i, L->x);
		ldl_dsolve(n, bp, D);
		ldl_ltsolve(n, bp, L->p, L->i, L->x);
		ldl_permt(n, x, bp, P);
	}
}
