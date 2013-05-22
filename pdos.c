#include "pdos.h"

// define "zero" threshold
static const double ZERO = 1e-8;

// constants and data structures
static const char* HEADER[] = {
  "Iter", 
  "||Ax+s-b||",
  " ||A'y+c||",
  "    c'x   ",
  "   -b'y   ",
  "    eta   ",
  "  lambda  "
};
// static const idxint LENS[] = {
//   4, 14, 11, 18, 11, 8, 9, 3
// }
static const idxint HEADER_LEN = 7;

// problem state
enum { SOLVED = 0, INDETERMINATE };

// to hold residual information
struct resid {
  double p_res;
  double d_res;
  double p_obj;
  double d_obj;
  double eta;
};

// forward declare inline declarations
static inline void relax(Work * w);
static inline void updateDualVars(Work * w);
static inline void adaptRhoAndSigma(Work * w, const idxint i);

static inline void prepZVariable(Work *w);
static inline void projectCones(Work * w,const Cone * k);
static inline void sety(const Work * w, Sol * sol);
static inline void sets(const Work * w, Sol * sol);
static inline void setx(const Work * w, Sol * sol);
static inline void getSolution(const Work * w, Sol * sol, idxint solver_state);
static inline void printSummary(const Work * w,idxint i, struct resid *r);
static inline void printHeader();
static inline void printSol(const Sol * sol);
static inline void freeWork(Work ** w);

Sol * pdos(const Data * d, const Cone * k)
{
  if(d == NULL || k == NULL) {
    return NULL;
  }
	idxint i, STATE = INDETERMINATE;
  struct resid residuals = { -1, -1, -1, -1, -1 };

  Params *p = d->p;
	Work * w = initWork(d, k);
  if(p->VERBOSE) {
    printHeader();
    tic();
  }
  
  for (i=0; i < p->MAX_ITERS; ++i){    
		projectLinSys(w);
		
    /* overrelaxation */
    relax(w);
    
    projectCones(w,k);
    updateDualVars(w);
    
    /* line search */
    // if(p->NORMALIZE) {
    //   adaptRhoAndSigma(w,i);
    // }
    
    residuals.p_res = calcPriResid(w);
    residuals.d_res = calcDualResid(w); 
    residuals.p_obj = calcPriObj(w);
    residuals.d_obj = calcDualObj(w);
    residuals.eta = fabs(residuals.p_obj - residuals.d_obj);
    
    
    // err = calcPriResid(d,w);
    // EPS_PRI = sqrt(w->l)*d->EPS_ABS + 
    //       d->EPS_REL*fmax(calcNorm(w->uv,w->l),calcNorm(w->uv_t,w->l));
    if (residuals.p_res < p->EPS_ABS && 
        residuals.d_res < p->EPS_ABS && 
        residuals.eta < p->EPS_ABS) {
      STATE = SOLVED;
      break;
    }
		if (p->VERBOSE && i % 10 == 0) printSummary(w,i, &residuals);
	}
	Sol * sol = PDOS_malloc(sizeof(Sol));
	getSolution(w,sol,STATE);
  
	if(p->VERBOSE) {
    printSummary(w,i,&residuals);
    PDOS_printf("Total solve time is %4.8fs\n", tocq());
	  //printSol(d,sol);
	}

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
  idxint i;
  for(i = 0; i < w->m; ++i) { w->y[i] += (w->s[i] - w->stilde[i])/w->lambda; }  
}

static inline void prepZVariable(Work *w){
  idxint i;
  for(i = 0; i < w->m; ++i) { w->s[i] = w->stilde[i] - w->lambda*w->y[i]; }
}

static inline void projectCones(Work * w,const Cone * k){
  // s = stilde - y
  // memcpy(w->s, w->stilde, d->m*sizeof(double));
  // addScaledArray(w->s, w->y, d->m, -1);
  // s = stilde - y
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
	//memcpy(sol->y, w->z + w->yi, d->m*sizeof(double));
  idxint i;
  for(i = 0; i < w->m; ++i) {
    sol->y[i] = w->D[i] * w->y[i];
  }
}

static inline void setx(const Work * w, Sol * sol){
  sol->n = w->n;
	sol->x = PDOS_malloc(sizeof(double)*w->n);
	//memcpy(sol->x, w->z, d->n*sizeof(double));
  idxint i;
  for(i = 0; i < w->n; ++i) {
    sol->x[i] = w->E[i] * w->x[i];
  }
}

static inline void sets(const Work * w, Sol * sol){
  sol->m = w->m;
	sol->s = PDOS_malloc(sizeof(double)*w->m);
	//memcpy(sol->y, w->z + w->yi, d->m*sizeof(double));
  idxint i;
  for(i = 0; i < w->m; ++i) {
    sol->s[i] = w->s[i] / w->D[i];
  }
}

static inline void relax(Work * w){   
	idxint j;
  const double ALPHA = w->params->ALPHA;
	for(j=0; j < w->m; ++j){
		w->stilde[j] = ALPHA*w->stilde[j] + (1.0 - ALPHA)*w->s[j];
	}  
}

static inline void adaptRhoAndSigma(Work *w, const idxint i) {
  const idxint MAX_ITERS = w->params->MAX_ITERS;
  const double TAU = w->params->TAU;
  const idxint SEARCH_ITERS = w->params->SEARCH_ITERS;
  const double BETA = w->params->BETA;
  
  if(i < MAX_ITERS/2) {
    // WARNING: major aliasing
    //   w->stilde is used for the temporary memory in both computations
    //   if we wanted to, we could allocate temp memory in the setup for this
    //   maybe it will speed it up... depends on compiler, architecture, et.c
    //
    
    // A'*y
    multByATrans(w,w->y,w->stilde);
    double norm_d = calcNormSq(w->stilde, w->n);
    if(norm_d < ZERO) {
      return;
    }
    double theta = -innerProd(w->stilde, w->c, w->n)/norm_d;
  
    // A*x + s
    memcpy(w->stilde, w->s, w->m*sizeof(double));
    accumByA(w,w->x,w->stilde);
    double norm_p = calcNormSq(w->stilde, w->m);
    if(norm_p < ZERO) {
      return;
    }
    double gamma = innerProd(w->stilde, w->b, w->m)/norm_p;
    
    if(theta > TAU && gamma > TAU && theta < 1.0/TAU && gamma < 1.0/TAU) {
      if (i <= SEARCH_ITERS) {
        w->lambda = 1.0/theta;
      } else {
        w->lambda = pow(w->lambda, 1.0-BETA)*pow(1.0/theta, BETA);
      }
    } else {
      return;
    }
  }
}

static inline void printSummary(const Work * w,idxint i, struct resid *r){
	// PDOS_printf("Iteration %i, primal residual %4f, primal tolerance %4f\n",i,err,EPS_PRI);
#ifdef DLONG
  PDOS_printf("%*li | ", (int)strlen(HEADER[0]), i);
#else
  PDOS_printf("%*i | ", (int)strlen(HEADER[0]), i);
#endif
  PDOS_printf("%*.3e   ", (int)strlen(HEADER[1]), r->p_res);
  PDOS_printf("%*.3e   ", (int)strlen(HEADER[2]), r->d_res);
  PDOS_printf("%*.3e   ", (int)strlen(HEADER[3]), r->p_obj);
  PDOS_printf("%*.3e   ", (int)strlen(HEADER[4]), r->d_obj);
  PDOS_printf("%*.3e   ", (int)strlen(HEADER[5]), r->eta);
  PDOS_printf("%*.3e\n", (int)strlen(HEADER[6]), w->lambda);  
}

static inline void printHeader() {
  idxint i, line_len;
  line_len = 0;
  for(i = 0; i < HEADER_LEN - 1; ++i) {
    PDOS_printf("%s | ", HEADER[i]);
    line_len += strlen(HEADER[i]) + 3;
  }
  PDOS_printf("%s\n", HEADER[HEADER_LEN-1]);
  line_len += strlen(HEADER[HEADER_LEN-1]);
  for(i = 0; i < line_len; ++i) {
    PDOS_printf("=");
  }
  PDOS_printf("\n");
}
