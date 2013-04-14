#include "private.h"
#include "common.h"

// forward declare
void choleskyInit(const cs * A, idxint P[], double **info);
void choleskyFactor(const cs * A, idxint P[], idxint Pinv[], cs ** L, double **D);
void choleskySolve(double *x, double b[], cs * L, double D[], idxint P[]);
void factorize(const Data * d,Work * w);

void freePriv(Work * w){
  cs_spfree(w->p->L);PDOS_free(w->p->P);PDOS_free(w->p->D);
  PDOS_free(w->p);
}

static inline void prepArgument(const Data *d, Work *w) {
  // memcpy(w->z_half,w->z,w->l*sizeof(double));
  // addScaledArray(w->z_half,w->lam,w->l,-1);
  
  // w->x = (-(x-c), b - (s+y))  
  idxint i;  
  for (i = 0; i < d->n; i++) { 
    // set x = -(x - c)
    w->x[i] = d->c[i]/w->rho - w->x[i];
  }
  for (i = 0; i < d->m; i++) { 
    // set s_half = (b - (s + y))
    w->stilde[i] = d->b[i]/w->sigma - (w->s[i] + w->y[i]);
  }
}

void projectLinSys(const Data *d, Work * w){
  // this preps the argument (-(x-c), b - (s+y)) for LDL solve 
  // (puts it in w->x)  
  prepArgument(d,w);  
  choleskySolve(w->x, w->x, w->p->L, w->p->D, w->p->P);  
  // stilde = b/sigma - A*x
  setAsScaledArray(w->stilde,d->b,(1/w->sigma),d->m); 
  // memcpy(w->stilde, d->b, d->m*sizeof(double));
  // stilde -= A*x
  decumByA(d, w->x, w->stilde);  
}

Work * initWork(const Data* d){ 
  idxint n_plus_m = d->n + d->m;
  Work *w = commonWorkInit(d);
  w->p = PDOS_malloc(sizeof(Priv));

  w->p->P = PDOS_malloc(sizeof(idxint)*n_plus_m);
  w->p->L = PDOS_malloc(sizeof (cs));
  w->p->L->m = n_plus_m;
  w->p->L->n = n_plus_m;
  w->p->L->nz = -1; 
  factorize(d,w);
  
  return w;
}

cs * formKKT(const Data * d, Work * w){
	/* ONLY UPPER TRIANGULAR PART IS STUFFED
	 * forms column compressed KKT matrix
	 * assumes column compressed form A matrix
   *
   * forms upper triangular part of [-I A'; A I]
	 */
	idxint j, k, kk;
	/* -I at top left */
  const idxint Anz = d->Ap[d->n];
	const idxint Knzmax = d->m + d->n + Anz;
	cs * K = cs_spalloc(d->m + d->n, d->m + d->n, Knzmax, 1, 1);
  kk = 0;
	for (k = 0; k < d->n; k++){
		K->i[kk] = k;
		K->p[kk] = k;
		K->x[kk] = -1;
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
	/* I at bottom right */
	for (k = 0; k < d->m; k++){
		K->i[kk] = k + d->n;
		K->p[kk] = k + d->n;
		K->x[kk] = 1;
    kk++;
	}
  // assert kk == Knzmax
	K->nz = Knzmax;
	cs * K_cs = cs_compress(K);
	cs_spfree(K);
  return(K_cs);
}


void factorize(const Data * d,Work * w){
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
  idxint * Pinv = cs_pinv(w->p->P, d->m+d->n);
  cs * C = cs_symperm(K, Pinv, 1); 
  choleskyFactor(C, NULL, NULL, &w->p->L, &w->p->D);
  if(d->VERBOSE) PDOS_printf("KKT matrix factorization took %4.8fs\n",tocq());
  cs_spfree(C);cs_spfree(K);PDOS_free(Pinv);PDOS_free(info);
}

void choleskyInit(const cs * A, idxint P[], double **info) {
	*info  = (double *) PDOS_malloc(AMD_INFO * sizeof(double));
#ifdef DLONG
	amd_l_order(A->n, A->p, A->i, P, (double *) NULL, *info);
#else
	amd_order(A->n, A->p, A->i, P, (double *) NULL, *info);
#endif
}

void choleskyFactor(const cs * A, idxint P[], idxint Pinv[], cs **L , double **D) 
{
	(*L)->p = (idxint *) PDOS_malloc((1 + A->n) * sizeof(idxint));
	idxint Parent[A->n], Lnz[A->n], Flag[A->n], Pattern[A->n];
	double Y[A->n];

#ifdef LDL_LONG
	ldl_l_symbolic(A->n, A->p, A->i, (*L)->p, Parent, Lnz, Flag, P, Pinv);
#else
	ldl_symbolic(A->n, A->p, A->i, (*L)->p, Parent, Lnz, Flag, P, Pinv);
#endif

	(*L)->nzmax = *((*L)->p + A->n);
  (*L)->x = (double *) PDOS_malloc((*L)->nzmax * sizeof(double));
	(*L)->i =    (idxint *) PDOS_malloc((*L)->nzmax * sizeof(idxint));
  *D  = (double *) PDOS_malloc(A->n * sizeof(double));
#ifdef LDL_LONG
	ldl_l_numeric(A->n, A->p, A->i, A->x, (*L)->p, Parent, Lnz, (*L)->i, (*L)->x, *D, Y, Pattern, Flag, P, Pinv);
#else
	ldl_numeric(A->n, A->p, A->i, A->x, (*L)->p, Parent, Lnz, (*L)->i, (*L)->x, *D, Y, Pattern, Flag, P, Pinv);
#endif
}

void choleskySolve(double *x, double b[], cs * L, double D[], idxint P[])
{
	if (P == NULL) {
		if (x != b) // if they're different addresses
      memcpy(x,b, (L->n)*sizeof(double)); 
    // for ( i = 0; i < n; i++) { x[i] = b[i]; } 

#ifdef LDL_LONG
    // do the solve
		ldl_l_lsolve(L->n, x, L->p, L->i, L->x);
		ldl_l_dsolve(L->n, x, D);
		ldl_l_ltsolve(L->n, x, L->p, L->i, L->x);
#else
    // do the solve
		ldl_lsolve(L->n, x, L->p, L->i, L->x);
		ldl_dsolve(L->n, x, D);
		ldl_ltsolve(L->n, x, L->p, L->i, L->x);
#endif  
	} else {
  	double bp1[L->n];

#ifdef LDL_LONG
    // do the solve
		ldl_l_perm(L->n, bp1, b, P);
		ldl_l_lsolve(L->n, bp1, L->p, L->i, L->x);
		ldl_l_dsolve(L->n, bp1, D);
		ldl_l_ltsolve(L->n, bp1, L->p, L->i, L->x);
		ldl_l_permt(L->n, x, bp1, P);
#else
    // do the solve
		ldl_perm(L->n, bp1, b, P);
		ldl_lsolve(L->n, bp1, L->p, L->i, L->x);
		ldl_dsolve(L->n, bp1, D);
		ldl_ltsolve(L->n, bp1, L->p, L->i, L->x);
		ldl_permt(L->n, x, bp1, P);
#endif
	}
}
