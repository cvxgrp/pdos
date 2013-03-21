#include "private.h"
#include "common.h"

Work * initWork(Data* d){
  
  int n_plus_m = d->n + d->m;
  Work * w = commonWorkInit(d);
  w->p = PDOS_malloc(sizeof(Priv));
  
  w->p->lambda = PDOS_calloc(n_plus_m + 1, sizeof(double));
  w->p->p = PDOS_calloc(n_plus_m + 1, sizeof(double));
  w->p->q = PDOS_calloc(n_plus_m + 1, sizeof(double));
  
  return w;
}

// void formQ(Data * d, Work * w){
//   int Qnzmax = 2*d->Ap[d->n]+2*d->n+2*d->m;
//   cs * Q = cs_spalloc(w->l/2, w->l/2, Qnzmax, 1, 1);
//   /* -A in Q : CCS: */
//   int j, k, kk=0;
//   for (j = 0; j < d->n; j++) {                 
//     for (k = d->Ap[j]; k < d->Ap[j+1]; k++) { 
//       Q->i[kk] = d->Ai[k] + d->n;
//       Q->p[kk] = j;
//       Q->x[kk] = -d->Ax[k];
//       kk++;
//     }   
//   }
//   /* A^T in Q : CCS: */
//   for (j = 0; j < d->n; j++) {                
//     for (k = d->Ap[j]; k < d->Ap[j+1]; k++) {
//       Q->p[kk] = d->Ai[k] + d->n;
//       Q->i[kk] = j;
//       Q->x[kk] = d->Ax[k];
//       kk++;
//     }   
//   }
//   /* -c^T in Q, same row */
//   for(k = 0; k < d->n; k++ ){
//     if (d->c[k]!=0){
//       Q->i[kk] = d->n + d->m;
//       Q->p[kk] = k;
//       Q->x[kk] = -d->c[k];
//       kk++;
//     }
//   }
//   /* c in Q, same col */
//   for(k = 0; k < d->n; k++ ){
//     if (d->c[k]!=0){
//       Q->i[kk] = k;
//       Q->p[kk] = d->n + d->m;
//       Q->x[kk] = d->c[k];
//       kk++;
//     }
//   }
//   /* -b^T in Q, same row */
//   for(k = 0; k < d->m; k++ ){
//     if (d->b[k]!=0){
//       Q->i[kk] = d->n + d->m;
//       Q->p[kk] = k + d->n;
//       Q->x[kk] = -d->b[k];
//       kk++;
//     }
//   }
//   /* b in Q, same col */
//   for(k = 0; k < d->m; k++ ){
//     if (d->b[k]!=0){
//       Q->i[kk] = k + d->n;
//       Q->p[kk] = d->n + d->m;
//       Q->x[kk] = d->b[k];
//       kk++;
//     }
//   }
//   Q->nz = kk;
//   w->p->Q = cs_compress(Q);
//   cs_spPDOS_free(Q);
// }

void freePriv(Work * w){
  PDOS_free(w->p->lambda); PDOS_free(w->p->p); PDOS_free(w->p->q); PDOS_free(w->p);
}

static inline void prepZHalfVariable(Data *d, Work *w) {
  // memcpy(w->z_half,w->z,w->l*sizeof(double));
  // addScaledArray(w->z_half,w->u,w->l,-1);
  
  // w->z_half = w->z - w->u = (x0,s0,r0,y0)  
  int i;  
  for (i = 0; i < d->n; i++) { 
    // set x_half = x - r_bar = x
    w->z_half[i] = w->z[i] - 0;
  }
  for (i = w->si; i < w->si + d->m; i++) { 
    // set s_half = (s - y_bar)
    w->z_half[i] =  w->z[i] - w->u[i];
  }
  for (i = w->ri; i < w->ri + d->n; i++) { 
    // set r_half = r - x_bar = -x_bar
    w->z_half[i] = 0 - w->u[i];
  }
  for (i = w->yi; i < w->yi + d->m; i++) {
    // set y_half = -(y - s_bar)
    w->z_half[i] = w->z[i] - w->u[i]; 
  }
}

// y += alpha*G*x
static inline void accumForward(const Data *d, const double *x, const double alpha, double *y) {
  // assumes x is an 2*(n+m) vector
  const double *input1 = x;
  const double *input2 = x + d->n;
  const double *input3 = x + d->n + d->m;
  const double *input4 = x + 2*(d->n) + d->m;
  
  // assumes y is an n+m+1 vector
  double *output1 = y;
  double *output2 = y + d->m;
  double *output3 = y + d->m + d->n; // a scalar
  
  /* output1 += alpha*A*input1 + alpha*input2 */
  // output1 += alpha*input2
  addScaledArray(output1, input2, d->m, alpha);
  // output1 += alpha*A*input1
  accumByScaledA(d, input1, alpha, output1);
  
  /* output2 += alpha*input3 - alpha*A'*input4 */
  // output2 += alpha*input3
  addScaledArray(output2, input3, d->n, alpha);
  // output2 += -alpha*A'*input4
  accumByScaledATrans(d, input4, -alpha, output2);
  
  /* output3 += alpha*(c'*input1 + b'*input4) */
  output3[0] += alpha*(innerProd(d->c, input1, d->n) + innerProd(d->b, input4, d->m));

}

// y += G^T*x
static inline void accumAdjoint(const Data *d, const double *x, double *y) {
  // assumes x is an n+m+1 vector
  const double *input1 = x;
  const double *input2 = x + d->m;
  const double input3 = x[d->m + d->n];
  
  // assumes y is an 2*(n+m) vector
  double *output1 = y;
  double *output2 = y + d->n;
  double *output3 = y + d->n + d->m;
  double *output4 = y + 2*(d->n) + d->m;
      
  /* output1 += A'*input1 + c*input3 */
  // output1 += c*input3
  addScaledArray(output1, d->c, d->n, input3);
  // output1 += A'*input1
  accumByATrans(d, input1, output1);
  
  /* output2 += input1 */
  addScaledArray(output2, input1, d->m, 1);

  /* output3 += input2 */
  addScaledArray(output3, input2, d->n, 1);
  
  /* output4 += -A*input2 + b*input3 */
  // output4 += b*input3
  addScaledArray(output4, d->b, d->m, input3);
  // output4 += -A*input2
  decumByA(d, input2, output4);
}

// y = G^T*x
static inline void applyAdjoint(const Data *d, const double *x, double *y) {
  // assumes x is an n+m+1 vector
  const double *input1 = x;
  const double *input2 = x + d->m;
  const double input3 = x[d->m + d->n];
  //PDOS_printf("%f\n",input3);
  
  // assumes y is an 2*(n+m) vector
  double *output1 = y;
  double *output2 = y + d->n;
  double *output3 = y + d->n + d->m;
  double *output4 = y + 2*(d->n) + d->m;
      
  /* output1 = A'*input1 + c*input3 */
  // output1 = c*input3
  setAsScaledArray(output1, d->c, input3, d->n);
  // output1 += A'*input1
  accumByATrans(d, input1, output1);
  
  /* output2 = input1 */
  memcpy(output2, input1, (d->m)*sizeof(double));
  // setAsScaledArray(output2, input1, 1, d->m);

  /* output3 = input2 */
  memcpy(output3, input2, (d->n)*sizeof(double));
  // setAsScaledArray(output3, input2, 1, d->n);
  
  /* output4 = -A*input2 + b*input3 */
  // output4 = b*input3
  setAsScaledArray(output4, d->b, input3, d->m);
  // output4 -= A*input2
  decumByA(d, input2, output4);
}


static inline void cgCustom(const Data *d, Work *w, int max_its,double tol){
	/* warm start cg with previous lambda */  
	int i = 0, n = d->n + d->m + 1;
  
  double *p = w->p->p;
  double *q = w->p->q;
  double *lambda = w->p->lambda;
  
	double alpha, beta, qsnew_sq=0;
  double tol_sq = tol*tol;  // XXX: could be a very small number...
  
  /* q = h - Gz_half - GG^T lambda */
  // q = h = (b,c,0)
  memcpy(q, d->b, (d->m)*sizeof(double));
  memcpy(q + d->m, d->c, (d->n)*sizeof(double));
  q[d->n + d->m] = 0;
  
  // for(i = 0; i < n; i++) {
  //   PDOS_printf("lamda[%d] = %f\n", i, lambda[i]);
  // }
  
  // w->ztmp = G^T*lambda (use previous lambda to warmstart)
  applyAdjoint(d, lambda, w->ztmp);
  
  // for(i = 0; i < w->l; i++) {
  //   PDOS_printf("G'lam[%d] = %f\n", i, w->ztmp[i]);
  // }
  
  // w->ztmp += w->z_half
  addScaledArray(w->ztmp, w->z_half, w->l, 1);
  // q += -G*w->ztmp
  accumForward(d, w->ztmp, -1, q);
  
  // for(i = 0; i < n; i++) {
  //   PDOS_printf("q[%d] = %f\n", i, q[i]);
  // }
  
  // p = q
	memcpy(p,q,n*sizeof(double));
  // ||q||^2
	double qsold_sq=calcNormSq(q,n);
  if (qsold_sq > tol_sq) {   // only iterate if the residual is small
  	for (i=0; i< max_its; ++i){
      // w->ztmp = G^T*p (multiplies by A and A^T)
      applyAdjoint(d, p, w->ztmp);
      // alpha = ||q||^2/ ||G^T*p||^2
  		alpha = qsold_sq/calcNormSq(w->ztmp, w->l);
    
      // lambda += alpha*p
      addScaledArray(lambda, p, n, alpha);
    
      // q += (-alpha)*G*w->ztmp
  		accumForward(d, w->ztmp, -alpha, q);    

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
	if (d->VERBOSE) PDOS_printf("terminating cg residual = %4f, took %i itns\n",sqrt(qsnew_sq),i);
}

void projectLinSys(Data * d, Work * w){
  // w->z_half = w->z - w->u
	memcpy(w->z_half,w->z,w->l*sizeof(double));
	addScaledArray(w->z_half,w->u,w->l,-1);    
  
  // w->ztmp = 0 (for temp computations)
  // memset(w->ztmp, 0, w->l*sizeof(double));

  // solve (GG^T)lambda = h - G*(w->z_half)
  cgCustom(d, w, d->CG_MAX_ITS, d->CG_TOL);
  
  // w->z_half += G'*lambda
  accumAdjoint(d, w->p->lambda, w->z_half);
}     


// void multByQ(const cs *Q, const double *x, double *y)
// {
//   /* y  = Q*x */
//   /* parallelizes over rows,
//      relies on special property that Q = -Q^T 
//      we use this so we can parallelize across columns
//      since Q is stored in compressed column format */
//   int p, j, n, *Qp, *Qi ;
//   double *Qx ;
//   n = Q->n ; Qp = Q->p ; Qi = Q->i ; Qx = Q->x ;
//   /* parallel matrix multiply */
//   int c1, c2;
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
