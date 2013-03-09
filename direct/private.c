#include "private.h"
// forward declare
void choleskyInit(cs * A, int P[], double **info);
void choleskyFactor(cs * A, int P[], int Pinv[], cs ** L, double **D);
void choleskySolve(double *x, double b[], cs * L, double D[], int P[]);
void factorize(Data * d,Work * w);

void freePriv(Work * w){
  cs_spfree(w->p->L);free(w->p->P);free(w->p->D);
  free(w->p->Ac); free(w->p->Atb);
  free(w->p->alpha1); free(w->p->alpha2);
  free(w->p);
}

static inline void prepZVariable(Data *d, Work *w) {
  // memcpy(w->z_half,w->z,w->l*sizeof(double));
  // addScaledArray(w->z_half,w->lam,w->l,-1);
  
  int i;
  // set x_half = x - r_bar = x
  for (i = 0; i < d->n; i++) { 
    w->z_half[i] = w->z[i] - 0; 
  }
  // set s_half = b - (s - y_bar)
  for (i = 0; i < d->m; i++) { 
    w->z_half[i + w->si] = d->b[i] - w->z[i + w->si] + w->u[i + w->si];
  }
  // set r_half = r - x_bar - c = -x_bar - c
  for (i = 0; i < d->n; i++) { 
    w->z_half[i + w->ri] = 0 - w->u[i + w->ri] - d->c[i];
  }
  // set y_half = -(y - s_bar)
  for (i = 0; i < d->n; i++) {
    w->z_half[i + w->yi] = - w->z[i + w->yi] + w->u[i + w->yi]; 
  }
}

void projectLinSys(Data *d, Work * w){

  prepZVariable(d,w);
  
  const double *x = w->z_half;
  const double *y = w->z_half + w->yi;

  choleskySolve(w->ztmp, w->z_half, w->p->L, w->p->D, w->p->P);
  
  double *kappa1 = (w->ztmp) + w->si; // has length m
  double *kappa2 = (w->ztmp) + w->ri; // has length n
  
  double delta3 = (innerProd(d->c,x,d->n) + innerProd(d->b,y,d->m) - \
    innerProd(w->p->Ac, kappa1, d->m) + innerProd(w->p->Atb, kappa2, d->n))/(w->p->s);
  
  // delta1 = kappa1 - delta3*alpha1
  // stored in (w->ztmp) + d->n
  addScaledArray(kappa1,w->p->alpha1, d->m,-delta3);
  // delta2 = kappa2 - delta3*alpha2
  // stored in (w->ztmp) + d->n + d->m
  addScaledArray(kappa2,w->p->alpha2, d->n,-delta3);
  
  // A'*delta1 + c*delta3
  setAsScaledArray(w->ztmp, d->c, delta3, d->n);
  accumByATrans(d, kappa1, w->ztmp);
  
  // -A*delta2 + b*delta3
  setAsScaledArray(w->ztmp + w->yi, d->b, delta3, d->m);
  decumByA(d, kappa1, w->ztmp);
  
  // z_half = z - [data.A'*delta1 + data.c*delta3; delta1; delta2; -data.A*delta2 + data.b*delta3]
  addScaledArray(w->z_half,w->ztmp,w->l,-1);
}

Work * initWork(Data* d){ 
  // CHECK
  int n_plus_m = d->n + d->m;
  Work * w = malloc(sizeof(Work));
  w->p = malloc(sizeof(Priv));
  w->l = 2*n_plus_m;
  w->si = d->n;
  w->ri = n_plus_m;
  w->yi = n_plus_m + d->n;
  w->z_half = malloc(sizeof(double)*w->l);
  w->z = calloc(w->l,sizeof(double));
  w->u = calloc(w->l,sizeof(double));  
  w->ztmp = calloc(w->l,sizeof(double));

  w->p->P = malloc(sizeof(int)*n_plus_m);
  w->p->L = malloc (sizeof (cs));
  w->p->L->m = n_plus_m;
  w->p->L->n = n_plus_m;
  w->p->L->nz = -1; 
  factorize(d,w);
  
  w->p->Ac = calloc(d->m, sizeof(double));
  w->p->Atb = calloc(d->n, sizeof(double));
  // Ac += A*c
  accumByA(d, d->c, w->p->Ac);
  // Atb += A'*b
  accumByATrans(d, d->b, w->p->Atb);
  
  double *tmp = calloc(w->l, sizeof(double));
  w->p->alpha1 = malloc(sizeof(double)*d->m);
  w->p->alpha2 = malloc(sizeof(double)*d->n);
  
  // memcpy(tmp, d->c, (d->n)*sizeof(double));
  setAsScaledArray(tmp, d->c, 1, d->n);
  // memcpy(tmp + 2*(d->n) + (d->m), d->b, (d->m)*sizeof(double));
  setAsScaledArray(tmp + 2*(d->n) + (d->m), d->b, -1, d->m);
  
  choleskySolve(w->ztmp, w->ztmp, w->p->L, w->p->D, w->p->P);
  
  // set alpha1 and alpha2
  memcpy(w->p->alpha1, tmp + (d->n), (d->m)*sizeof(double));
  memcpy(w->p->alpha2, tmp + (d->n) + (d->m), (d->n)*sizeof(double));
  
  free(tmp);

  w->p->s = calcNormSq(d->c, d->n) + calcNormSq(d->b, d->m) - \
    innerProd(w->p->Ac, w->p->alpha1, d->m) + innerProd(w->p->Atb, w->p->alpha2, d->n);

  return w;
}

cs * formKKT(Data * d, Work * w){
	/* ONLY UPPER TRIANGULAR PART IS STUFFED
	 * forms column compressed KKT matrix
	 * assumes column compressed form A matrix
   *
   * forms upper triangular part of [I A'; A -I]
	 */
	int j, k;
	/* I at top left */
  const int Anz = d->Ap[d->n];
	const int Knzmax = w->l + 2*Anz;
	cs * K = cs_spalloc(d->m + d->n, d->m + d->n, Knzmax, 1, 1);
	for (k = 0; k < d->n; k++){
		K->i[k] = k;
		K->p[k] = k;
		K->x[k] = 1;
	}
	/* A^T at top right : CCS: */
	for (j = 0; j < d->n; j++) {                 
		for (k = d->Ap[j]; k < d->Ap[j+1]; k++) { 
			K->p[k] = d->Ai[k] + d->n;
			K->i[k] = j;
			K->x[k] = d->Ax[k];
		}   
	}
  
	/* -I at bottom right */
	for (k = 0; k < d->m; k++){
		K->i[k] = k + d->n;
		K->p[k] = k + d->n;
		K->x[k] = -1;
	}
	K->nz = Knzmax;
	cs * K_cs = cs_compress(K);
	cs_spfree(K);
  return(K_cs);
}


void factorize(Data * d,Work * w){
  tic();
  cs * K = formKKT(d,w);
  printf("KKT matrix factorization info:\n");
  double *info;
  choleskyInit(K, w->p->P, &info);
  amd_info(info);
  int * Pinv = cs_pinv(w->p->P, w->l);
  cs * C = cs_symperm(K, Pinv, 1); 
  choleskyFactor(C, NULL, NULL, &w->p->L, &w->p->D);
  printf("KKT matrix factorization took %4.2f s\n",tocq());
  cs_spfree(C);cs_spfree(K);free(Pinv);free(info);
}

void choleskyInit(cs * A, int P[], double **info) {
	*info  = (double *) malloc(AMD_INFO * sizeof(double));
	amd_order(A->n, A->p, A->i, P, (double *) NULL, *info);
}

void choleskyFactor(cs * A, int P[], int Pinv[], cs **L , double **D) 
{
	(*L)->p = (int *) malloc((1 + A->n) * sizeof(int));
	int Parent[A->n], Lnz[A->n], Flag[A->n], Pattern[A->n];
	double Y[A->n];

	ldl_symbolic(A->n, A->p, A->i, (*L)->p, Parent, Lnz, Flag, P, Pinv);

	(*L)->nzmax = *((*L)->p + A->n);
  (*L)->x = (double *) malloc((*L)->nzmax * sizeof(double));
	(*L)->i =    (int *) malloc((*L)->nzmax * sizeof(int));
  *D  = (double *) malloc(A->n * sizeof(double));

	ldl_numeric(A->n, A->p, A->i, A->x, (*L)->p, Parent, Lnz, (*L)->i, (*L)->x, *D, Y, Pattern, Flag, P, Pinv);
}

void choleskySolve(double *x, double b[], cs * L, double D[], int P[])
{
  // solves two systems at once (use openMP here?)
  // have to trust that x and b have length 2n
  // solves Ax1 = b1, Ax2 = b2 for x = (x1,x2), b = (b1,b2)
  // int i;
	if (P == NULL) {
		if (x != b) // if they're different addresses
      memcpy(x,b, 2*(L->n)*sizeof(double)); 
    // for ( i = 0; i < n; i++) { x[i] = b[i]; } 
    
    // do the first solve   
		ldl_lsolve(L->n, x, L->p, L->i, L->x);
		ldl_dsolve(L->n, x, D);
		ldl_ltsolve(L->n, x, L->p, L->i, L->x);
    
    // do the second solve
		ldl_lsolve(L->n, x + (L->n), L->p, L->i, L->x);
		ldl_dsolve(L->n, x + (L->n), D);
		ldl_ltsolve(L->n, x + (L->n), L->p, L->i, L->x);
	} else {
  	double bp1[L->n];
    double bp2[L->n];
    
    // do the first solve
		ldl_perm(L->n, bp1, b, P);
		ldl_lsolve(L->n, bp1, L->p, L->i, L->x);
		ldl_dsolve(L->n, bp1, D);
		ldl_ltsolve(L->n, bp1, L->p, L->i, L->x);
		ldl_permt(L->n, x, bp1, P);
    
    // do the second solve
		ldl_perm(L->n, bp2, b + (L->n), P);
		ldl_lsolve(L->n, bp2, L->p, L->i, L->x);
		ldl_dsolve(L->n, bp2, D);
		ldl_ltsolve(L->n, bp2, L->p, L->i, L->x);
		ldl_permt(L->n, x + (L->n), bp2, P);
	}
}
