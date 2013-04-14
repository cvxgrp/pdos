#include "pdos.h"

// define "zero" threshold
static const double ZERO = 1e-8;

// constants and data structures
static const char* HEADER[] = {
  "Iter", 
  "ni(Ax+s-b)",
  " ni(A'y+c)",
  "    c'u   ",
  "   -b'w   ",
  "    eta   ",
  "    rho   ",
  "   sigma  "
};
// static const idxint LENS[] = {
//   4, 14, 11, 18, 11, 8, 9, 3
// }
static const idxint HEADER_LEN = 8;

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
static inline void relax(const Data * d, Work * w);
static inline void updateDualVars(const Data *d, Work * w);
static inline void adaptRhoAndSigma(const Data * d, Work * w, const idxint i);

static inline void prepZVariable(const Data *d, Work *w);
static inline void projectCones(const Data * d,Work * w,const Cone * k);
static inline void sety(const Data * d, const Work * w, Sol * sol);
static inline void setx(const Data * d, const Work * w, Sol * sol);
static inline void getSolution(const Data * d, const Work * w, Sol * sol, idxint solver_state);
static inline void printSummary(const Data * d,const Work * w,idxint i, struct resid *r);
static inline void printHeader();
static inline void printSol(const Data * d, const Sol * sol);
static inline void freeWork(Work * w);

Sol * pdos(const Data * d, const Cone * k)
{
  if(d == NULL || k == NULL) {
    return NULL;
  }
	idxint i, STATE = INDETERMINATE;
  struct resid residuals = { -1, -1, -1, -1, -1 };

	Work * w = initWork(d);
  if(d->VERBOSE) {
    printHeader();
    tic();
  }
  for (i=0; i < d->MAX_ITERS; ++i){    
		projectLinSys(d,w);
		
    /* overrelaxation */
    relax(d,w);
    
    projectCones(d,w,k);
    updateDualVars(d,w);
    
    /* line search */
    if(d->NORMALIZE) {
      adaptRhoAndSigma(d,w,i);
    }
    
    residuals.p_res = calcPriResid(d,w);
    residuals.d_res = calcDualResid(d,w); 
    residuals.p_obj = calcPriObj(d,w);
    residuals.d_obj = calcDualObj(d,w);
    residuals.eta = fabs(residuals.p_obj - residuals.d_obj);
    
    
    // err = calcPriResid(d,w);
    // EPS_PRI = sqrt(w->l)*d->EPS_ABS + 
    //       d->EPS_REL*fmax(calcNorm(w->uv,w->l),calcNorm(w->uv_t,w->l));
    if (residuals.p_res < d->EPS_ABS*w->dual_scale && 
        residuals.d_res < d->EPS_ABS*w->primal_scale && 
        residuals.eta < d->EPS_ABS) {
      STATE = SOLVED;
      break;
    }
		if (d->VERBOSE && i % 10 == 0) printSummary(d,w,i, &residuals);
	}
	Sol * sol = PDOS_malloc(sizeof(Sol));
	getSolution(d,w,sol,STATE);
  
	if(d->VERBOSE) {
    printSummary(d,w,i,&residuals);
    PDOS_printf("Total solve time is %4.8fs\n", tocq());
	  //printSol(d,sol);
	}
  if(d->NORMALIZE) {
    // unscale data
    for(i = 0; i < d->Ap[d->n]; ++i) {
      d->Ax[i] /= (w->dual_scale)*(w->primal_scale);
    }
    for(i = 0; i < d->n; ++i) {
      d->c[i] /= w->primal_scale;
    }
    for(i = 0; i < d->m; ++i) {
      d->b[i] /= w->dual_scale;
    }
  }
  freeWork(w);
	return sol;
}

void free_data(Data * d, Cone * k){
  if(d) {
    if(d->b) PDOS_free(d->b);
    if(d->c) PDOS_free(d->c);
    if(d->Ax) PDOS_free(d->Ax);
    if(d->Ai) PDOS_free(d->Ai);
    if(d->Ap) PDOS_free(d->Ap);
    PDOS_free(d);
  }
  if(k) {
    if(k->q) PDOS_free(k->q);
    PDOS_free(k);
  }
  d = NULL; k = NULL;
}

void free_sol(Sol *sol){
  if(sol) {
    if(sol->x) PDOS_free(sol->x);
    if(sol->y) PDOS_free(sol->y);
    // done automatically
    // if(sol->status) PDOS_free(sol->status);
    PDOS_free(sol);
  }
  sol = NULL;
}

static inline void freeWork(Work * w){
  freePriv(w);
  PDOS_free(w->x); // also frees w->stilde
  w->stilde = NULL;
  PDOS_free(w->s);
  PDOS_free(w->y);
  PDOS_free(w);
}

static inline void printSol(const Data * d, const Sol * sol){
	idxint i;
	PDOS_printf("%s\n",sol->status); 
	if (sol->x != NULL){
		for ( i=0;i<d->n; ++i){
#ifdef DLONG
			PDOS_printf("x[%li] = %4f\n",i, sol->x[i]);
#else
			PDOS_printf("x[%i] = %4f\n",i, sol->x[i]);
#endif
		}
	}
	if (sol->y != NULL){
		for ( i=0;i<d->m; ++i){
#ifdef DLONG
			PDOS_printf("y[%li] = %4f\n",i, sol->y[i]);
#else
			PDOS_printf("y[%i] = %4f\n",i, sol->y[i]);
#endif
		}
	}
}

static inline void updateDualVars(const Data *d, Work * w){
  idxint i;
  for(i = 0; i < d->m; ++i) { w->y[i] += (w->s[i] - w->stilde[i]); }  
}

static inline void prepZVariable(const Data *d, Work *w){
  idxint i;
  for(i = 0; i < d->m; ++i) { w->s[i] = w->stilde[i] - w->y[i]; }
}

static inline void projectCones(const Data *d,Work * w,const Cone * k){
  // s = stilde - y
  // memcpy(w->s, w->stilde, d->m*sizeof(double));
  // addScaledArray(w->s, w->y, d->m, -1);
  // s = stilde - y
  prepZVariable(d,w);

	/* s onto K */
	projCone(w->s, k);
}

static inline void getSolution(const Data * d, const Work * w, Sol * sol, idxint solver_state){
  setx(d,w,sol);
  sety(d,w,sol);
  switch(solver_state) {
    case SOLVED: memcpy(sol->status,"Solved", 7*sizeof(char)); break;
    default: memcpy(sol->status, "Indeterminate", 15*sizeof(char));
  }
}

static inline void sety(const Data * d,const Work * w, Sol * sol){
	sol->y = PDOS_malloc(sizeof(double)*d->m);
	//memcpy(sol->y, w->z + w->yi, d->m*sizeof(double));
  idxint i;
  for(i = 0; i < d->m; ++i) {
    sol->y[i] = w->rho * w->dual_scale * w->y[i];
  }
}

static inline void setx(const Data * d,const Work * w, Sol * sol){
	sol->x = PDOS_malloc(sizeof(double)*d->n);
	//memcpy(sol->x, w->z, d->n*sizeof(double));
  idxint i;
  for(i = 0; i < d->n; ++i) {
    sol->x[i] = w->sigma * w->primal_scale * w->x[i];
  }
}

static inline void relax(const Data * d,Work * w){   
	idxint j;  
	for(j=0; j < d->m; ++j){
		w->stilde[j] = (d->ALPH)*w->stilde[j] + (1.0 - d->ALPH)*w->s[j];
	}  
}

static inline void adaptRhoAndSigma(const Data *d, Work *w, const idxint i) {
  if(i < d->MAX_ITERS/2) {
    // WARNING: major aliasing
    //   w->stilde is used for the temporary memory in both computations
    //   if we wanted to, we could allocate temp memory in the setup for this
    //   maybe it will speed it up... depends on compiler, architecture, et.c
    //
    
    // A'*y
    multByATrans(d,w->y,w->stilde);
    double norm_d = calcNormSq(w->stilde, d->n);
    if(norm_d < ZERO) {
      return;
    }
    double theta = -innerProd(w->stilde, d->c, d->n)/norm_d;
  
    // A*x + s
    memcpy(w->stilde, w->s, d->m*sizeof(double));
    accumByA(d,w->x,w->stilde);
    double norm_p = calcNormSq(w->stilde, d->m);
    if(norm_p < ZERO) {
      return;
    }
    double gamma = innerProd(w->stilde, d->b, d->m)/norm_p;
    
    if(theta > d->TAU && gamma > d->TAU && theta < 1.0/d->TAU && gamma < 1.0/d->TAU) {
      if (i <= d->SEARCH_ITERS) {
        w->rho = theta;
        w->sigma = gamma;
      } else {
        w->rho = pow(w->rho, 1.0-d->BETA)*pow(theta, d->BETA);
        w->sigma = pow(w->sigma,1.0-d->BETA)*pow(gamma, d->BETA);
      }
    } else {
      return;
    }
  }
}

static inline void printSummary(const Data * d,const Work * w,idxint i, struct resid *r){
	// PDOS_printf("Iteration %i, primal residual %4f, primal tolerance %4f\n",i,err,EPS_PRI);
#ifdef DLONG
  PDOS_printf("%*li | ", (int)strlen(HEADER[0]), i);
#else
  PDOS_printf("%*i | ", (int)strlen(HEADER[0]), i);
#endif
  PDOS_printf("%*.4f   ", (int)strlen(HEADER[1]), r->p_res);
  PDOS_printf("%*.4f   ", (int)strlen(HEADER[2]), r->d_res);
  PDOS_printf("%*.4f   ", (int)strlen(HEADER[3]), r->p_obj);
  PDOS_printf("%*.4f   ", (int)strlen(HEADER[4]), r->d_obj);
  PDOS_printf("%*.4f   ", (int)strlen(HEADER[5]), r->eta);
  PDOS_printf("%*.4f   ", (int)strlen(HEADER[5]), w->rho);
  PDOS_printf("%*.4f\n", (int)strlen(HEADER[5]), w->sigma);
  
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
