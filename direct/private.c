#include "private.h"
#include "common.h"

/* forward declare */
void choleskyInit(const cs * A, idxint P[], double **info);
void choleskyFactor(const cs * A, idxint P[], idxint Pinv[], cs ** L, double **D);
void choleskySolve(double *x, double b[], cs * L, double D[], idxint P[]);
void factorize(Work * w);

void freePriv(Work * w){
  cs_spfree(w->p->L);PDOS_free(w->p->P);PDOS_free(w->p->D);
  PDOS_free(w->p);
}

static __inline void prepArgument(Work *w) {
  /* uses the x and stilde memory space to store the argument, since we don't
   * need that memory space anymore
   */

  /* w->x = (-(x-lambda*c), b - (s+y)) */
  idxint i;
  for (i = 0; i < w->n; i++) {
    /* set x = -(x - lambda*c) */
    w->x[i] = w->lambda*w->c[i] - w->x[i];
  }
  for (i = 0; i < w->m; i++) {
    /* set stilde = (b - (s + lambda*y)) */
    w->stilde[i] = w->b[i] - (w->s[i] + w->lambda*w->y[i]);
  }
}

void projectLinSys(Work * w){
  idxint i = 0;

  /* this preps the argument (-(x-lambda*c), b - (s+lambda*y)) for LDL solve
   * (puts it in w->x)
   */
  prepArgument(w);

  /* now solve the linear system */
  choleskySolve(w->x, w->x, w->p->L, w->p->D, w->p->P);

  /* now find stilde... */
  /* stilde = b - A*x */
  /* now, stilde contains b - v - A*x,
  * where v = s + lambda * y. 
  * To recover stilde = b - A*x, we simply add "v"
  * to the solution.
  */
  /* stilde = b - A*x */
  for (i = 0; i < w->m; i++) {
   /* set stilde = (b + (s + lambda*y)) */
   w->stilde[i] += (w->s[i] + w->lambda*w->y[i]);
  }
}

Work * initWork(const Data* d, const Cone *k){
  /* initialize workspace for direct solver */
  idxint n_plus_m = d->n + d->m;
  Work *w = commonWorkInit(d,k);
  w->p = PDOS_malloc(sizeof(Priv));

  /* allocate permutation vector and LDL data structures */
  w->p->P = PDOS_malloc(sizeof(idxint)*n_plus_m);
  w->p->L = PDOS_malloc(sizeof (cs));
  w->p->L->m = n_plus_m;
  w->p->L->n = n_plus_m;
  w->p->L->nz = -1;

  /* factorize the KKT system */
  factorize(w);

  return w;
}

cs * formKKT(Work * w){
  /* ONLY UPPER TRIANGULAR PART IS STUFFED
   * forms column compressed KKT matrix
   * assumes column compressed form A matrix
   *
   * forms upper triangular part of [-I A'; A I]
   */
  idxint j, k, kk;
  /* -I at top left */
  const idxint Anz = w->Ap[w->n];
  const idxint Knzmax = w->m + w->n + Anz;
  cs * K = cs_spalloc(w->m + w->n, w->m + w->n, Knzmax, 1, 1);
  cs * K_cs = NULL;

  kk = 0;
  for (k = 0; k < w->n; k++){
    K->i[kk] = k;
    K->p[kk] = k;
    K->x[kk] = -1;
    kk++;
  }

  /* A^T at top right : CCS: */
  for (j = 0; j < w->n; j++) {
    for (k = w->Ap[j]; k < w->Ap[j+1]; k++) {
      K->p[kk] = w->Ai[k] + w->n;
      K->i[kk] = j;
      K->x[kk] = w->Ax[k];
      kk++;
    }
  }
  /* I at bottom right */
  for (k = 0; k < w->m; k++){
    K->i[kk] = k + w->n;
    K->p[kk] = k + w->n;
    K->x[kk] = 1;
    kk++;
  }
  /* assert kk == Knzmax */
  K->nz = Knzmax;
  K_cs = cs_compress(K);
  cs_spfree(K);
  return(K_cs);
}

void factorize(Work * w){
  cs * K = NULL;
  static timer KKT_timer;
  double *info;
  idxint * Pinv;
  cs * C;

  tic(&KKT_timer);

  K = formKKT(w);
#ifdef PRINTKKT
  if(w->params->VERBOSE) PDOS_printf("Factorization info:\n");
#endif
  choleskyInit(K, w->p->P, &info);

  /* perform ordering */
#ifdef PRINTKKT
  if(w->params->VERBOSE) {
#ifdef DLONG
    amd_l_info(info);
#else
    amd_info(info);
#endif
  }
#endif

  /* compute the inverse permutation */
  Pinv = cs_pinv(w->p->P, w->m+w->n);

  /* permute the KKT matrix */
  C = cs_symperm(K, Pinv, 1);

  /* perform the LDL factorization */
  choleskyFactor(C, NULL, NULL, &w->p->L, &w->p->D);

  if(w->params->VERBOSE) PDOS_printf("Factorization time: %4.8fs\n",tocq(&KKT_timer));
  cs_spfree(C);cs_spfree(K);
  PDOS_free(Pinv);PDOS_free(info);
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
  idxint *Parent = (idxint *) PDOS_malloc(A->n * sizeof(idxint));
  idxint *Lnz = (idxint *) PDOS_malloc(A->n * sizeof(idxint));
  idxint *Flag = (idxint *) PDOS_malloc(A->n * sizeof(idxint));
  idxint *Pattern = (idxint *) PDOS_malloc(A->n * sizeof(idxint));
  double *Y = (double *) PDOS_malloc(A->n * sizeof(double));

  (*L)->p = (idxint *) PDOS_malloc((1 + A->n) * sizeof(idxint));

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
  PDOS_free(Parent); PDOS_free(Lnz); PDOS_free(Flag);
  PDOS_free(Pattern); PDOS_free(Y);
}

void choleskySolve(double *x, double b[], cs * L, double D[], idxint P[])
{
  /* return x = P^{-T} L^{-T} D^{-1} L^{-1} P^{-1} b */
	if (P == NULL) {
    /* if x = b, then will do an in-place solve */
		if (x != b) /* if they're different addresses */
      memcpy(x,b, (L->n)*sizeof(double));

#ifdef LDL_LONG
    /* do the solve */
		ldl_l_lsolve(L->n, x, L->p, L->i, L->x);
		ldl_l_dsolve(L->n, x, D);
		ldl_l_ltsolve(L->n, x, L->p, L->i, L->x);
#else
    /* do the solve */
		ldl_lsolve(L->n, x, L->p, L->i, L->x);
		ldl_dsolve(L->n, x, D);
		ldl_ltsolve(L->n, x, L->p, L->i, L->x);
#endif
	} else {
  	double *bp1 = malloc(L->n * sizeof(double));

#ifdef LDL_LONG
    /* do the solve */
		ldl_l_perm(L->n, bp1, b, P);
		ldl_l_lsolve(L->n, bp1, L->p, L->i, L->x);
		ldl_l_dsolve(L->n, bp1, D);
		ldl_l_ltsolve(L->n, bp1, L->p, L->i, L->x);
		ldl_l_permt(L->n, x, bp1, P);
#else
    /* do the solve */
		ldl_perm(L->n, bp1, b, P);
		ldl_lsolve(L->n, bp1, L->p, L->i, L->x);
		ldl_dsolve(L->n, bp1, D);
		ldl_ltsolve(L->n, bp1, L->p, L->i, L->x);
		ldl_permt(L->n, x, bp1, P);
#endif
    free(bp1);
	}
}
