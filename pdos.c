#include "pdos.h"

// redefine printfs and memory allocators as needed
#ifdef MATLAB_MEX_FILE
#include "mex.h"
int (*pdos_printf) (const char *, ...) = mexPrintf ;
void (*pdos_free) (void *) = mxFree ;
void *(*pdos_malloc) (size_t) = mxMalloc ;
void *(*pdos_calloc) (size_t, size_t) = mxCalloc ;
#elif defined PYTHON
#include <Python.h>
#include <stdlib.h>
int (*pdos_printf) (const char *, ...) = PySys_WriteStdout ;
void (*pdos_free) (void *) = free ;
void *(*pdos_malloc) (size_t) = malloc ;
void *(*pdos_calloc) (size_t, size_t) = calloc ;
#else
#include <stdio.h>
#include <stdlib.h>
int (*pdos_printf) (const char *, ...) = printf ;
void (*pdos_free) (void *) = free ;
void *(*pdos_malloc) (size_t) = malloc ;
void *(*pdos_calloc) (size_t, size_t) = calloc ;
#endif

// constants and data structures
static const char* HEADER[] = {
  "Iter", 
  "ni(Ax+s-b)",
  " ni(A'y+c)",
  " ni(Au+v) ",
  " ni(A'w)  ",
  "    c'u   ",
  "   -b'w   ",
  "    eta   "
};
// static const int LENS[] = {
//   4, 14, 11, 18, 11, 8, 9, 3
// }
static const int HEADER_LEN = 8;

// problem state
enum { SOLVED = 0, INFEASIBLE, UNBOUNDED, INDETERMINATE };

// to hold residual information
struct resid {
  double p_res;
  double d_res;
  double p_inf;
  double d_inf;
  double p_obj;
  double d_obj;
  double eta;
};

// forward declare inline declarations
static inline void relax(Data * d, Work * w);
static inline void updateDualVars(Work * w);
static inline void prepZVariable(Work *w);
static inline void projectCones(Data * d,Work * w,Cone * k);
static inline void sety(Data * d, Work * w, Sol * sol);
static inline void setx(Data * d, Work * w, Sol * sol);
static inline void getSolution(Data * d, Work * w, Sol * sol, int solver_state);
static inline void printSummary(Data * d,Work * w,int i, struct resid *r);
static inline void printHeader();
static inline void printSol(Data * d, Sol * sol);
static inline void freeWork(Work * w);

Sol * pdos(Data * d, Cone * k)
{
  if(d == NULL || k == NULL) {
    return NULL;
  }
	int i, STATE = INDETERMINATE;
  struct resid residuals = { -1, -1, -1, -1, -1, -1, -1 };

	Work * w = initWork(d);
  if(d->VERBOSE) {
    printHeader();
    tic();
  }
  for (i=0; i < d->MAX_ITERS; ++i){             
		projectLinSys(d,w);
		relax(d,w);
    projectCones(d,w,k);
    updateDualVars(w);
    
    residuals.p_res = calcPriResid(d,w);
    residuals.d_res = calcDualResid(d,w);
    residuals.eta = calcSurrogateGap(d,w);
        
    residuals.p_inf = calcCertPriResid(d,w);
    residuals.d_inf = calcCertDualResid(d,w);
    residuals.p_obj = calcCertPriObj(d,w);
    residuals.d_obj = calcCertDualObj(d,w);
    
    // err = calcPriResid(d,w);
    // EPS_PRI = sqrt(w->l)*d->EPS_ABS + 
    //       d->EPS_REL*fmax(calcNorm(w->uv,w->l),calcNorm(w->uv_t,w->l));
    if (residuals.p_res < d->EPS_ABS*w->dual_scale && 
        residuals.d_res < d->EPS_ABS*w->primal_scale && 
        residuals.eta < d->EPS_ABS) {
      STATE = SOLVED;
      break;
    }
    if (residuals.d_inf < d->EPS_INFEAS*w->dual_scale && 
        residuals.p_inf < d->EPS_INFEAS*w->primal_scale && 
        residuals.p_obj - residuals.d_obj < -d->EPS_INFEAS /* fix at 1e-2? relative? */) {
      if (residuals.d_obj < d->EPS_INFEAS)
        STATE = UNBOUNDED;
      else if (-residuals.p_obj < d->EPS_INFEAS)
        STATE = INFEASIBLE;
      else
        STATE = INDETERMINATE;
      break;
    }
		if (d->VERBOSE && i % 10 == 0) printSummary(d,w,i, &residuals);
	}
	Sol * sol = pdos_malloc(sizeof(Sol));
	getSolution(d,w,sol,STATE);
  
	if(d->VERBOSE) {
    printSummary(d,w,i,&residuals);
    pdos_printf("Total solve time is %4.8fs\n", tocq());
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
    if(d->b) pdos_free(d->b);
    if(d->c) pdos_free(d->c);
    if(d->Ax) pdos_free(d->Ax);
    if(d->Ai) pdos_free(d->Ai);
    if(d->Ap) pdos_free(d->Ap);
    pdos_free(d);
  }
  if(k) {
    if(k->q) pdos_free(k->q);
    pdos_free(k);
  }
  d = NULL; k = NULL;
}

void free_sol(Sol *sol){
  if(sol) {
    if(sol->x) pdos_free(sol->x);
    if(sol->y) pdos_free(sol->y);
    if(sol->status) pdos_free(sol->status);
    pdos_free(sol);
  }
  sol = NULL;
}

static inline void freeWork(Work * w){
  freePriv(w);
  pdos_free(w->z_half);
  pdos_free(w->z);
  pdos_free(w->u);
  pdos_free(w->ztmp);
  pdos_free(w);
}

static inline void printSol(Data * d, Sol * sol){
	int i;
	pdos_printf("%s\n",sol->status); 
	if (sol->x != NULL){
		for ( i=0;i<d->n; ++i){
			pdos_printf("x[%i] = %4f\n",i, sol->x[i]);
		}
	}
	if (sol->y != NULL){
		for ( i=0;i<d->m; ++i){
			pdos_printf("y[%i] = %4f\n",i, sol->y[i]);
		}
	}
}

static inline void updateDualVars(Work * w){
  int i;
  for(i = 0; i < w->l; ++i) { w->u[i] += (w->ztmp[i] - w->z[i]); }
}

static inline void prepZVariable(Work *w){
  int i;
  for(i = 0; i < w->l; ++i) { w->z[i] = w->ztmp[i] + w->u[i]; }
}

static inline void projectCones(Data *d,Work * w,Cone * k){
  // memcpy(w->z, w->ztmp, w->l*sizeof(double));
  // addScaledArray(w->z, w->u, w->l,1);
  
  // z = z_half + u
  prepZVariable(w);
  
  /* x onto R^n */
  // do nothing
  
	/* s onto K */
	projCone(w->z + w->si, k);

	/* r onto 0 */
	memset(w->z + w->ri, 0, (sizeof(double)*d->n));
  
  /* y onto K^* */
	projDualCone(w->z + w->yi, k);
}

static inline void getSolution(Data * d, Work * w, Sol * sol, int solver_state){
  setx(d,w,sol);
  sety(d,w,sol);
  switch(solver_state) {
    case SOLVED: sol->status = strdup("Solved"); break;
    case INFEASIBLE: sol->status = strdup("Infeasible"); break;
    case UNBOUNDED: sol->status = strdup("Unbounded"); break;
    default: sol->status = strdup("Indeterminate");
  }
}

static inline void sety(Data * d,Work * w, Sol * sol){
	sol->y = pdos_malloc(sizeof(double)*d->m);
	//memcpy(sol->y, w->z + w->yi, d->m*sizeof(double));
  int i;
  for(i = 0; i < d->m; ++i) {
    sol->y[i] = w->dual_scale * w->z[i + w->yi];
  }
}

static inline void setx(Data * d,Work * w, Sol * sol){
	sol->x = pdos_malloc(sizeof(double)*d->n);
	//memcpy(sol->x, w->z, d->n*sizeof(double));
  int i;
  for(i = 0; i < d->n; ++i) {
    sol->x[i] = w->primal_scale * w->z[i];
  }
}

static inline void relax(Data * d,Work * w){   
	int j;
	for(j=0; j < w->l; ++j){
		w->ztmp[j] = d->ALPH*w->z_half[j] + (1.0 - d->ALPH)*w->z[j];
	}
}

static inline void printSummary(Data * d,Work * w,int i, struct resid *r){
	// pdos_printf("Iteration %i, primal residual %4f, primal tolerance %4f\n",i,err,EPS_PRI);
  pdos_printf("%*i | ", (int)strlen(HEADER[0]), i);
  pdos_printf("%*.4f   ", (int)strlen(HEADER[1]), r->p_res); // p_res
  pdos_printf("%*.4f   ", (int)strlen(HEADER[2]), r->d_res); // d_res
  pdos_printf("%*.4f   ", (int)strlen(HEADER[3]), r->p_inf); // full(p_inf));
  pdos_printf("%*.4f   ", (int)strlen(HEADER[4]), r->d_inf);//full(d_inf));
  pdos_printf("%*.4f   ", (int)strlen(HEADER[5]), r->p_obj);//full(p_obj));
  pdos_printf("%*.4f   ", (int)strlen(HEADER[6]), r->d_obj);//full(d_obj));
  pdos_printf("%*.4f\n", (int)strlen(HEADER[7]), r->eta);//full(eta));
}

static inline void printHeader() {
  int i, line_len;
  line_len = 0;
  for(i = 0; i < HEADER_LEN - 1; ++i) {
    pdos_printf("%s | ", HEADER[i]);
    line_len += strlen(HEADER[i]) + 3;
  }
  pdos_printf("%s\n", HEADER[HEADER_LEN-1]);
  line_len += strlen(HEADER[HEADER_LEN-1]);
  for(i = 0; i < line_len; ++i) {
    pdos_printf("=");
  }
  pdos_printf("\n");
}
