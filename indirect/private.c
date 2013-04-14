#include "private.h"
#include "common.h"

Work * initWork(const Data* d){
  
  Work * w = commonWorkInit(d);
  w->p = PDOS_malloc(sizeof(Priv));
  
  w->p->Ax = PDOS_calloc(d->m, sizeof(double));
  w->p->p = PDOS_calloc(d->n, sizeof(double));
  w->p->q = PDOS_calloc(d->n, sizeof(double));
  
  return w;
}

void freePriv(Work * w){
  PDOS_free(w->p->p); PDOS_free(w->p->q); PDOS_free(w->p->Ax); PDOS_free(w->p);
}

static inline void prepArgument(const Data *d, Work *w) {
  // memcpy(w->z_half,w->z,w->l*sizeof(double));
  // addScaledArray(w->z_half,w->lam,w->l,-1);
  
  // w->x = (x, s + y - b)
  
  idxint i;  
  // for (i = 0; i < d->n; i++) { 
  //   // set x = (x - c)
  //   w->x[i] -= d->c[i]/w->rho;
  // }
  for (i = 0; i < d->m; i++) { 
    // set s_half = (s + y - b/sigma)
    w->stilde[i] = w->s[i] + w->y[i] - d->b[i]/w->sigma;
  }
}

static inline void cgCustom(const Data *d, Work *w, idxint max_its,double tol){
	/* warm start cg with previous x */  
	idxint i = 0;
  const idxint n = d->n;
  
  double *p = w->p->p;
  double *q = w->p->q;
  double *Ax = w->p->Ax;
  double *x = w->x; // contains x
  const double *s = w->stilde; // contains v - b/sigma
  
	double alpha, beta, qsnew_sq=0;
  double tol_sq = tol*tol;  // XXX: could be a very small number...
  
  /* q = -c - A'*(A*x - b + v) */
  memcpy(Ax, s, (d->m)*sizeof(double));
  setAsScaledArray(q,d->c,-(1.0/w->rho),d->n);  // q = -c/rho
  accumByA(d,x,Ax);   // Ax = A*x + (v - b/sigma)
  decumByATrans(d,Ax,q); // q = -c/rho - A'*(A*x - b/sigma + v)
  
  // p = q
	memcpy(p,q,n*sizeof(double));
  // ||q||^2
	double qsold_sq=calcNormSq(q,n);
  if (qsold_sq > tol_sq) {   // only iterate if the residual is small
  	for (i=0; i< max_its; ++i){
      // uses stilde as temp variable
      // stilde = p + A'*A*p
      multByA(d,p,Ax); // Ax = A*p
      memcpy(w->stilde, p, n*sizeof(double));
      accumByATrans(d,Ax,w->stilde);
      
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
	if (d->VERBOSE) PDOS_printf("terminating cg residual = %4f, took %li itns\n",sqrt(qsnew_sq),i);
#else
	if (d->VERBOSE) PDOS_printf("terminating cg residual = %4f, took %i itns\n",sqrt(qsnew_sq),i);
#endif
}

void projectLinSys(const Data * d, Work * w){
  // this only modifies w->s = s + y - b
  prepArgument(d,w);

  // solve (I+A'*A)*x = u + A'*(b - v)
  cgCustom(d, w, d->CG_MAX_ITS, d->CG_TOL);
  
  // stilde = b/sigma - A*x
  setAsScaledArray(w->stilde,d->b,(1.0/w->sigma),d->m); 
  
  // stilde -= A*x
  decumByA(d, w->x, w->stilde);
}     


// void multByQ(const cs *Q, const double *x, double *y)
// {
//   /* y  = Q*x */
//   /* parallelizes over rows,
//      relies on special property that Q = -Q^T 
//      we use this so we can parallelize across columns
//      since Q is stored in compressed column format */
//   idxint p, j, n, *Qp, *Qi ;
//   double *Qx ;
//   n = Q->n ; Qp = Q->p ; Qi = Q->i ; Qx = Q->x ;
//   /* parallel matrix multiply */
//   idxint c1, c2;
//   double yj;
// #pragma omp parallel for private(p,c1,c2,yj) 
//   for (j = 0 ; j < n ; j++)
//   {
//     yj = 0.0;
//     c1 = Qp[j]; c2 = Qp[j+1];
//     for (p = c1 ; p < c2 ; p++)        
//     {   
//       yj -= Qx[p] * x[ Qi[p] ] ;
//     }
//     y[j] = yj;
//   }
// }
