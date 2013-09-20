#include "pdos.h"
#ifndef NDEBUG
#include <assert.h>
#endif

// define "zero" threshold
static const double ZERO = 1e-8;

// constants and data structures
static const char* HEADER[] = {
  " iter",
  " ||Ax+s-b||",
  " ||A'y+c||",
  "    c'x    ",
  "   -b'y    ",
  "    eta   "
};
// static const idxint LENS[] = {
//   4, 14, 11, 18, 11, 8, 9, 3
// }
static const idxint HEADER_LEN = 6;

// problem state
enum { SOLVED = 0, INDETERMINATE };

// to hold residual information
struct resid {
  double p_res;
  double d_res;
  double p_obj;
  double d_obj;
  double eta;
  double eps_gap;
};

// forward declare inline declarations
static inline void relax(Work * w);
static inline void updateDualVars(Work * w);
//static inline void adaptRhoAndSigma(Work * w, const idxint i);

static inline void prepZVariable(Work *w);
static inline void projectCones(Work * w,const Cone * k);
static inline void sety(const Work * w, Sol * sol);
static inline void sets(const Work * w, Sol * sol);
static inline void setx(const Work * w, Sol * sol);
static inline void getSolution(const Work * w, Sol * sol, idxint solver_state);
static inline void printSummary(idxint i, struct resid *r);
static inline void printHeader();
static inline void printSol(const Sol * sol);
static inline void freeWork(Work ** w);

Sol * pdos(const Data * d, const Cone * k)
{
  static timer PDOS_timer;

  if(d == NULL || k == NULL) {
    return NULL;
  }
  idxint i, STATE = INDETERMINATE;
  struct resid residuals = { -1, -1, -1, -1, -1 };

#ifndef NDEBUG
  // ensure that cone sizes match data size
  idxint cone_sz = 0;
  for (i = 0; i < k->qsize; ++i) {
    cone_sz += k->q[i];
  }
  assert( d->m == cone_sz + k->f + k->l );
#endif

  // set the parameters
  Params *p = d->p;
  // set the denominators of DIMACS error measures (i.e., relative scale of
  // error)
  const double pscale = (1.0 + calcNormInf(d->b, d->m));
  const double dscale = (1.0 + calcNormInf(d->c, d->n));

  if(p->VERBOSE) {
    PDOS_printf("\nPDOS - A Primal-Dual Operator Splitting for Cone Programming.\n");
    PDOS_printf("       (c) E. Chu, B. O'Donoghue, N. Parikh, S. Boyd, Stanford University, 2012-13.\n\n");
  }

  // initialize workspace, allocates memory for necessary computations
  Work * w = initWork(d, k);
  if(!w->p) {
    PDOS_printf("Error with Pardiso license.\n");
    freeWork(&w);
    return NULL;
  }

  if(p->VERBOSE) {
    PDOS_printf("lambda: %5.3e\n\n", w->lambda);

    printHeader();
    tic(&PDOS_timer);
  }

  for (i=0; i < p->MAX_ITERS; ++i){
    // Pi_P
	projectLinSys(w);

    /* overrelaxation */
    relax(w);

    // Pi_K
    projectCones(w,k);
    // y += (1.0/lambda)*(s - stilde)
    updateDualVars(w);

    residuals.p_res = calcPriResid(w) ;
    residuals.d_res = calcDualResid(w);
    residuals.p_obj = calcPriObj(w);
    residuals.d_obj = calcDualObj(w);
    residuals.eta = fabs(residuals.p_obj - residuals.d_obj);

    // check against DIMACS error measures
    if (residuals.p_res < p->EPS_ABS * pscale &&
        residuals.d_res < p->EPS_ABS * dscale &&
        residuals.eta < p->EPS_ABS * (1.0 + fabs(residuals.p_obj) + fabs(residuals.d_obj))) {
      STATE = SOLVED;
      break;
    }
		if (p->VERBOSE && i % 10 == 0) printSummary(i, &residuals);
	}
	Sol * sol = PDOS_malloc(sizeof(Sol));
	getSolution(w,sol,STATE);

	if(p->VERBOSE) {
      printSummary(i,&residuals);
      PDOS_printf("Total solve time is %4.8fs\n", tocq(&PDOS_timer));
	}

  // free the temporary workspace
  freeWork(&w);
	return sol;
}

void freeData(Data **d, Cone **k){
  if(*d) {
    if((*d)->b) PDOS_free((*d)->b);
    if((*d)->c) PDOS_free((*d)->c);
    if((*d)->Ax) PDOS_free((*d)->Ax);
    if((*d)->Ai) PDOS_free((*d)->Ai);
    if((*d)->Ap) PDOS_free((*d)->Ap);
    if((*d)->p) PDOS_free((*d)->p);
    PDOS_free(*d);
  }
  if(*k) {
    if((*k)->q) PDOS_free((*k)->q);
    PDOS_free(*k);
  }
  *d = NULL; *k = NULL;
}

void freeSol(Sol **sol){
  if(*sol) {
    if((*sol)->x) PDOS_free((*sol)->x);
    if((*sol)->y) PDOS_free((*sol)->y);
    if((*sol)->s) PDOS_free((*sol)->s);
    // done automatically
    // if(sol->status) PDOS_free(sol->status);
    PDOS_free(*sol);
  }
  *sol = NULL;
}

static inline void freeWork(Work **w){
  if(*w) {
    freePriv(*w);
    if((*w)->x) PDOS_free((*w)->x); // also frees w->stilde
    (*w)->stilde = NULL;
    if((*w)->s) PDOS_free((*w)->s);
    if((*w)->y) PDOS_free((*w)->y);
    if((*w)->params->NORMALIZE) { // if normalized, we have new memory
      if((*w)->Ax) PDOS_free((*w)->Ax);
      if((*w)->b) PDOS_free((*w)->b);
      if((*w)->c) PDOS_free((*w)->c);
    }
    if((*w)->Atx) PDOS_free((*w)->Atx);
    if((*w)->Ati) PDOS_free((*w)->Ati);
    if((*w)->Atp) PDOS_free((*w)->Atp);

    if((*w)->D) PDOS_free((*w)->D);
    if((*w)->E) PDOS_free((*w)->E);
    PDOS_free(*w);
  }
  *w = NULL;
}

static inline void printSol(const Sol * sol){
	idxint i;
	PDOS_printf("%s\n",sol->status);
	if (sol->x != NULL){
		for ( i=0;i< sol->n; ++i){
#ifdef DLONG
			PDOS_printf("x[%li] = %4f\n",i, sol->x[i]);
#else
			PDOS_printf("x[%i] = %4f\n",i, sol->x[i]);
#endif
		}
	}
	if (sol->y != NULL){
		for ( i=0;i<sol->m; ++i){
#ifdef DLONG
			PDOS_printf("y[%li] = %4f\n",i, sol->y[i]);
#else
			PDOS_printf("y[%i] = %4f\n",i, sol->y[i]);
#endif
		}
	}
	if (sol->s != NULL){
		for ( i=0;i<sol->m; ++i){
#ifdef DLONG
			PDOS_printf("s[%li] = %4f\n",i, sol->s[i]);
#else
			PDOS_printf("s[%i] = %4f\n",i, sol->s[i]);
#endif
		}
	}
}

static inline void updateDualVars(Work * w){
  // y = y + (1/lambda)*(s - stilde)
  idxint i;
  for(i = 0; i < w->m; ++i) { w->y[i] += (w->s[i] - w->stilde[i])/w->lambda; }
}

static inline void prepZVariable(Work *w){
  idxint i;
  for(i = 0; i < w->m; ++i) { w->s[i] = w->stilde[i] - w->lambda*w->y[i]; }
}

static inline void projectCones(Work * w,const Cone * k){
  // s = stilde - lambda*y
  prepZVariable(w);

	/* s onto K */
	projCone(w->s, k);
}

static inline void getSolution(const Work * w, Sol * sol, idxint solver_state){
  setx(w,sol);
  sety(w,sol);
  sets(w,sol);
  switch(solver_state) {
    case SOLVED: memcpy(sol->status,"Solved", 7*sizeof(char)); break;
    default: memcpy(sol->status, "Indeterminate", 15*sizeof(char));
  }
}

static inline void sety(const Work * w, Sol * sol){
  sol->m = w->m;
	sol->y = PDOS_malloc(sizeof(double)*w->m);

  // y = D*y (scale the variable)
  idxint i;
  for(i = 0; i < w->m; ++i) {
    sol->y[i] = w->D[i] * w->y[i];
  }
}

static inline void setx(const Work * w, Sol * sol){
  sol->n = w->n;
	sol->x = PDOS_malloc(sizeof(double)*w->n);

  // x = E*x (scale the variable)
  idxint i;
  for(i = 0; i < w->n; ++i) {
    sol->x[i] = w->E[i] * w->x[i];
  }
}

static inline void sets(const Work * w, Sol * sol){
  sol->m = w->m;
	sol->s = PDOS_malloc(sizeof(double)*w->m);

  // s = D^{-1}*s (scale the variable)
  idxint i;
  for(i = 0; i < w->m; ++i) {
    sol->s[i] = w->s[i] / w->D[i];
  }
}

static inline void relax(Work * w){
  // stilde = alpha*stilde + (1 - alpha)*s
	idxint j;
  const double ALPHA = w->params->ALPHA;
	for(j=0; j < w->m; ++j){
		w->stilde[j] = ALPHA*w->stilde[j] + (1.0 - ALPHA)*w->s[j];
	}
}


static inline void printSummary(idxint i, struct resid *r){
#ifdef DLONG
  PDOS_printf("%*li | ", (int)strlen(HEADER[0]), i);
#else
  PDOS_printf("%*i | ", (int)strlen(HEADER[0]), i);
#endif
  PDOS_printf("%*.3e   ", (int)strlen(HEADER[1]), r->p_res);
  PDOS_printf("%*.3e   ", (int)strlen(HEADER[2]), r->d_res);
  PDOS_printf("%*.3e   ", (int)strlen(HEADER[3]), r->p_obj);
  PDOS_printf("%*.3e   ", (int)strlen(HEADER[4]), r->d_obj);
  PDOS_printf("%*.3e   \n", (int)strlen(HEADER[5]), r->eta);
}

static inline void printHeader() {
  idxint i, line_len;
  line_len = 0;
  for(i = 0; i < HEADER_LEN - 1; ++i) {
    if (i == 0 )
      PDOS_printf("%s | ", HEADER[i]);
    else
      PDOS_printf("%s   ", HEADER[i]);
    line_len += strlen(HEADER[i]) + 3;
  }
  PDOS_printf("%s\n", HEADER[HEADER_LEN-1]);
  line_len += strlen(HEADER[HEADER_LEN-1]);
  for(i = 0; i < line_len+3; ++i) {
    PDOS_printf("=");
  }
  PDOS_printf("\n");
}
