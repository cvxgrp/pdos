#include "private.h"
#include "common.h"

// forward declare
void choleskyInit(cs * A, int P[], double **info);
void choleskyFactor(cs * A, int P[], int Pinv[], cs ** L, double **D);
void choleskySolve(double *x, double b[], cs * L, double D[], int P[]);
void factorize(Data * d,Work * w);

void freePriv(Work * w){
  cs_spfree(w->p->L);PDOS_free(w->p->P);PDOS_free(w->p->D);
  PDOS_free(w->p->Ac); PDOS_free(w->p->Atb);
  PDOS_free(w->p->alpha1); PDOS_free(w->p->alpha2);
  PDOS_free(w->p);
}

static inline void prepZHalfVariable(Data *d, Work *w) {
  // memcpy(w->z_half,w->z,w->l*sizeof(double));
  // addScaledArray(w->z_half,w->lam,w->l,-1);
  
  // w->z_half = w->z - w->u = (x0,s0,r0,y0)
  // w->ztmp, on the other hand, will contain (x0,b-s0,-c+r0,-y0)
  
  int i;  
  for (i = 0; i < d->n; i++) { 
    // set x_half = x - r_bar = x
    w->z_half[i] = w->z[i] - 0;
    // set x = x - r_bar = x
    w->ztmp[i] = w->z_half[i];
  }
  for (i = w->si; i < w->si + d->m; i++) { 
    // set s_half = (s - y_bar)
    w->z_half[i] =  w->z[i] - w->u[i];
    // set s = b - (s - y_bar)
    w->ztmp[i] = d->b[i - w->si] - w->z_half[i];
  }
  for (i = w->ri; i < w->ri + d->n; i++) { 
    // set r_half = r - x_bar = -x_bar
    w->z_half[i] = 0 - w->u[i];
    // set r = r - x_bar - c = -x_bar - c
    w->ztmp[i] = w->z_half[i] - d->c[i - w->ri];
  }
  for (i = w->yi; i < w->yi + d->m; i++) {
    // set y_half = -(y - s_bar)
    w->z_half[i] = w->z[i] - w->u[i]; 
    // set y = -(y - s_bar)
    w->ztmp[i] = -w->z_half[i]; 
  }
}

void projectLinSys(Data *d, Work * w){
  // this preps z_half = z - u, and ztmp = argument for LDL solve
  prepZHalfVariable(d,w);
  
  const double *x = w->z_half;
  const double *y = w->z_half + w->yi;  
  
  choleskySolve(w->ztmp, w->ztmp, w->p->L, w->p->D, w->p->P);
  
  double *kappa1 = (w->ztmp) + w->si; // has length m
  double *kappa2 = (w->ztmp) + w->ri; // has length n
  
  // (c'*x + b'*y - (A*c)'*kappa1 + (A'*b)'*kappa2)/s
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
  decumByA(d, kappa2, w->ztmp + w->yi);
  
  // z_half = z - [data.A'*delta1 + data.c*delta3; delta1; delta2; -data.A*delta2 + data.b*delta3]
  addScaledArray(w->z_half,w->ztmp,w->l,-1);
}

Work * initWork(Data* d){ 
  int n_plus_m = d->n + d->m;
  Work *w = commonWorkInit(d);
  w->p = PDOS_malloc(sizeof(Priv));

  w->p->P = PDOS_malloc(sizeof(int)*n_plus_m);
  w->p->L = PDOS_malloc(sizeof (cs));
  w->p->L->m = n_plus_m;
  w->p->L->n = n_plus_m;
  w->p->L->nz = -1; 
  factorize(d,w);
  
  w->p->Ac = PDOS_calloc(d->m, sizeof(double));
  w->p->Atb = PDOS_calloc(d->n, sizeof(double));
  // Ac += A*c
  accumByA(d, d->c, w->p->Ac);
  // Atb += A'*b
  accumByATrans(d, d->b, w->p->Atb);
  
  // tmp may be extraneous (put into w->ztmp?)
  double *tmp = PDOS_calloc(w->l, sizeof(double));
  w->p->alpha1 = PDOS_malloc(sizeof(double)*d->m);
  w->p->alpha2 = PDOS_malloc(sizeof(double)*d->n);
  
  // memcpy(tmp, d->c, (d->n)*sizeof(double));
  setAsScaledArray(tmp, d->c, 1, d->n);
  // memcpy(tmp + 2*(d->n) + (d->m), d->b, (d->m)*sizeof(double));
  setAsScaledArray(tmp + 2*(d->n) + (d->m), d->b, -1, d->m);
  
  choleskySolve(tmp, tmp, w->p->L, w->p->D, w->p->P);
  
  // set alpha1 and alpha2
  memcpy(w->p->alpha1, tmp + (d->n), (d->m)*sizeof(double));
  memcpy(w->p->alpha2, tmp + (d->n) + (d->m), (d->n)*sizeof(double));
  
  PDOS_free(tmp);

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
	int j, k, kk;
	/* I at top left */
  const int Anz = d->Ap[d->n];
	const int Knzmax = (w->l)/2 + Anz;
	cs * K = cs_spalloc(d->m + d->n, d->m + d->n, Knzmax, 1, 1);
  kk = 0;
	for (k = 0; k < d->n; k++){
		K->i[kk] = k;
		K->p[kk] = k;
		K->x[kk] = 1;
    kk++;
	}
	/* A^T at top right : CCS: */
	for (j = 0; j < d->n; j++) {                 
		for (k = d->Ap[j]; k < d->Ap[j+1]; k++) { 
			K->p[kk] = d->Ai[k] + d->n;
			K->i[kk] = j;
			K->x[kk] = d->Ax[k];
      kk++;
		}   
	}
	/* -I at bottom right */
	for (k = 0; k < d->m; k++){
		K->i[kk] = k + d->n;
		K->p[kk] = k + d->n;
		K->x[kk] = -1;
    kk++;
	}
  // assert kk == Knzmax
	K->nz = Knzmax;
	cs * K_cs = cs_compress(K);
	cs_spfree(K);
  return(K_cs);
}


void factorize(Data * d,Work * w){
  tic();
  cs * K = formKKT(d,w);
  if(d->VERBOSE) PDOS_printf("KKT matrix factorization info:\n");
  double *info;
  choleskyInit(K, w->p->P, &info);
  if(d->VERBOSE) {
#ifdef DLONG
    amd_l_info(info);
#else
    amd_info(info);
#endif
  }
  int * Pinv = cs_pinv(w->p->P, (w->l)/2);
  cs * C = cs_symperm(K, Pinv, 1); 
  choleskyFactor(C, NULL, NULL, &w->p->L, &w->p->D);
  if(d->VERBOSE) PDOS_printf("KKT matrix factorization took %4.8fs\n",tocq());
  cs_spfree(C);cs_spfree(K);PDOS_free(Pinv);PDOS_free(info);
}

void choleskyInit(cs * A, int P[], double **info) {
	*info  = (double *) PDOS_malloc(AMD_INFO * sizeof(double));
#ifdef DLONG
	amd_l_order(A->n, A->p, A->i, P, (double *) NULL, *info);
#else
	amd_order(A->n, A->p, A->i, P, (double *) NULL, *info);
#endif
}

void choleskyFactor(cs * A, int P[], int Pinv[], cs **L , double **D) 
{
	(*L)->p = (int *) PDOS_malloc((1 + A->n) * sizeof(int));
	int Parent[A->n], Lnz[A->n], Flag[A->n], Pattern[A->n];
	double Y[A->n];

	ldl_symbolic(A->n, A->p, A->i, (*L)->p, Parent, Lnz, Flag, P, Pinv);

	(*L)->nzmax = *((*L)->p + A->n);
  (*L)->x = (double *) PDOS_malloc((*L)->nzmax * sizeof(double));
	(*L)->i =    (int *) PDOS_malloc((*L)->nzmax * sizeof(int));
  *D  = (double *) PDOS_malloc(A->n * sizeof(double));

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
