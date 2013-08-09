#include "private.h"
#include "common.h"

Work * initWork(const Data *d, const Cone *k){
  // private indirect method initialization
  Work * w = commonWorkInit(d,k);
  w->p = PDOS_malloc(sizeof(Priv));

  // space for vectors containing A*x and for CG's p and q
  w->p->Ax = PDOS_calloc(d->m, sizeof(double));
  w->p->p = PDOS_calloc(d->n, sizeof(double));
  w->p->q = PDOS_calloc(d->n, sizeof(double));

  return w;
}

void freePriv(Work * w){
  PDOS_free(w->p->p); PDOS_free(w->p->q); PDOS_free(w->p->Ax); PDOS_free(w->p);
}

static inline void prepArgument(Work *w) {
  // uses the stilde memory space to store the argument, since we don't need
  // that memory space anymore
  idxint i;
  for (i = 0; i < w->m; i++) {
    // set stilde = s + y - b
    w->stilde[i] = w->s[i] + w->lambda*w->y[i] - w->b[i];
  }
}

static inline void cgCustom(Work *w){
  // solve (I+A^TA) x = (x^k - lambda c + A^T*(b - s^k - lambda y^k))
  // recall that s^k = b - A*x^k
  const idxint MAX_ITERS = w->params->CG_MAX_ITS;
  const double TOL = w->params->CG_TOL;
  const idxint VERBOSE = w->params->VERBOSE;

	/* warm start cg with previous x */
	idxint i = 0;
  const idxint n = w->n;

  double *p = w->p->p;
  double *q = w->p->q;
  double *Ax = w->p->Ax;
  double *x = w->x; // contains x
  const double *s = w->stilde; // contains v - b

	double alpha, beta, qsnew_sq=0;
  // we multiply the tolerance by RATIO since the "x" space is scaled by the
  // inverse of this ratio
  double tol_sq = TOL*TOL*RATIO;  // XXX: could be a very small number...

  /* q = -lambda * c - A'*(A*x - b + v) */
  memcpy(Ax, s, (w->m)*sizeof(double));
  setAsScaledArray(q,w->c,-w->lambda,w->n);  // q = -c*lambda
  accumByA(w,x,Ax);   // Ax = A*x + (v - b)
  decumByATrans(w,Ax,q); // q = -c*lambda - A'*(A*x - b + v)

  // p = q
	memcpy(p,q,n*sizeof(double));
  // ||q||^2
	double qsold_sq=calcNormSq(q,n);
  if (qsold_sq > tol_sq) {   // only iterate if the residual is small
  	for (i=0; i< MAX_ITERS; ++i){
      // uses stilde as temp variable
      // stilde = p + A'*A*p
      multByA(w,p,Ax); // Ax = A*p
      memcpy(w->stilde, p, n*sizeof(double));
      accumByATrans(w,Ax,w->stilde);

      // alpha = ||q||^2/ (p'*(I + A'*A)*p)
  		alpha = qsold_sq/innerProd(p,w->stilde,n);

      // x += alpha*p
      addScaledArray(x, p, n, alpha);

      // q += (-alpha)*(p + A'*A*p)
      addScaledArray(q, w->stilde, n, -alpha);

      // ||q||^2
  		qsnew_sq = calcNormSq(q,n);
  		if (qsnew_sq < tol_sq){
        ++i;
  			break;
  		}

      // beta = ||q_{i+1}||^2 / ||q_i||^2
      beta = qsnew_sq / qsold_sq;
      // p = q + beta*p
  		scaleArray(p,beta,n);     // p = beta*p
      addScaledArray(p,q,n,1);  // p += q
  		qsold_sq = qsnew_sq;
  	}
  }
#ifdef DLONG
  if (VERBOSE) PDOS_printf("terminating cg residual = %4f, took %li itns\n",sqrt(qsnew_sq),i);
#else
  if (VERBOSE) PDOS_printf("terminating cg residual = %4f, took %i itns\n",sqrt(qsnew_sq),i);
#endif
}

void projectLinSys(Work * w){
  // this only modifies w->s = s + y - b
  prepArgument(w);

  // solve (I+A'*A)*x = u + A'*(b - v)
  cgCustom(w);

  // stilde = b - A*x
  memcpy(w->stilde, w->b, w->m*sizeof(double));

  // stilde -= A*x
  decumByA(w, w->x, w->stilde);
}

