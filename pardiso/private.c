#ifndef NDEBUG
#include <assert.h>
#endif

#include "private.h"
#include "common.h"

/* PARDISO prototype. */
void pardisoinit (void   *, int *,   int *, int *, double *, int *);
void pardiso     (void   *, int *,   int *, int *,    int *, int *,
                  double *, int *,   int *, int *,   int *, int *,
                  int *, double *, double *, int *, double *);
void pardiso_chkmatrix    (int *, int *, double *, int *, int *, int *);
void pardiso_chkvec       (int *, int *, double *, int *);
void pardiso_printstats   (int *, int *, double *, int *, int *, int *,
                           double *, int *);


void formKKT(Work *w);

Work * initWork(const Data *d, const Cone *k){
  static timer init_timer;

  // private indirect method initialization
  tic(&init_timer);
  Work * w = commonWorkInit(d,k);
  PDOS_printf("Common work initialization took %0.4fs\n", tocq(&init_timer));

  tic(&init_timer);
  w->p = PDOS_malloc(sizeof(Priv));

  /* set up Pardiso parameters */
  w->p->mtype = -2;
  w->p->error = 0;
  w->p->solver = 0; /* use sparse direct solver */
  // w->p->solver = 1; /* use iterative solver */
  w->p->n = w->m + w->n;
  w->p->tmp = PDOS_calloc(w->p->n, sizeof(double));

  /* Numbers of processors, value of OMP_NUM_THREADS */
  int num_procs;
  char *var = getenv("OMP_NUM_THREADS");
  if(var != NULL)
    sscanf( var, "%d", &num_procs );
  else {
    PDOS_printf("Set environment OMP_NUM_THREADS to 1\n");
    exit(1);
  }

  w->p->iparm[2] = num_procs; /* number of processors */  

  w->p->maxfct = 1;		/* Maximum number of numerical factorizations.  */
  w->p->mnum   = 1;         /* Which factorization to use. */
    
  w->p->msglvl = 1;         /* Print statistical information  */

  formKKT(w); 

  pardisoinit (w->p->pt,  &(w->p->mtype), &(w->p->solver), w->p->iparm, w->p->dparm, &(w->p->error)); 

  // w->p->dparm[0] = 20; // MAX CG ITS
  
  PDOS_printf("Pardiso init took %0.4fs\n", tocq(&init_timer));

  if (w->p->error != 0) 
  {
    if (w->p->error == -10 )
      PDOS_printf("No license file found \n");
    if (w->p->error == -11 )
      PDOS_printf("License is expired \n");
    if (w->p->error == -12 )
      PDOS_printf("Wrong username or hostname \n");
    freePriv(w);
    return w; 
  }
  else
    PDOS_printf("[PARDISO]: License check was successful ... \n");

  int nrhs = 1;
  int idum; double ddum; 
#ifndef NDEBUG 
  tic(&init_timer);
  double *b = PDOS_calloc(w->p->n, sizeof(double));
  pardiso_printstats(&(w->p->mtype), &(w->p->n), w->p->U, w->p->Ui, w->p->Uj, &nrhs, b, &(w->p->error));
  PDOS_printf("Checking stats took %0.4fs\n", tocq(&init_timer));
  PDOS_printf("  ERROR CODE: %d\n", w->p->error);
  PDOS_free(b);
#endif

  /* now, perform numerical factorization */
  tic(&init_timer); 
/* -------------------------------------------------------------------- */
/* ..  Reordering and Numeric Factorization.  This step also allocates  */
/*     all memory that is necessary for the factorization.              */
/* -------------------------------------------------------------------- */
  w->p->phase = 12; 
  //w->p->phase = 11; // if iterative 
  w->p->iparm[27] = 1;  /* parallel metis reordering if 1 */   
  pardiso (w->p->pt, &(w->p->maxfct), &(w->p->mnum), &(w->p->mtype), 
           &(w->p->phase), &(w->p->n), w->p->U, w->p->Ui, w->p->Uj, 
           &idum, &nrhs, w->p->iparm, &(w->p->msglvl), &ddum, &ddum, 
           &(w->p->error), w->p->dparm);

  if (w->p->error != 0) {
      PDOS_printf("ERROR during symbolic factorization: %d\n", w->p->error);
      freePriv(w);
      return w;
  }
  PDOS_printf("Reordering completed ... \n");
  PDOS_printf("Number of nonzeros in factors  = %d\n", w->p->iparm[17]);
  PDOS_printf("Number of factorization MFLOPS = %d\n", w->p->iparm[18]);
  PDOS_printf("Factorization time: %0.4fs\n\n", tocq(&init_timer));

  w->p->msglvl = 0; // make quiet

  // get ready to solve
  w->p->phase = 33;
  
  w->p->iparm[5] = 1;       /* write solution to b (rhs) */
  w->p->iparm[7] = 0;       /* Max numbers of iterative refinement steps. */
  w->p->iparm[20] = 0;      /* only use 1x1 diagonal pivoting */

  return w;
}

void freePriv(Work * w) {
  int idum;
  double ddum;
  int nrhs = 1;
/* -------------------------------------------------------------------- */    
/* ..  Termination and release of memory.                               */
/* -------------------------------------------------------------------- */    
  w->p->phase = -1;                 /* Release internal memory. */
    
  
  if (w->p->error == 0) {
    pardiso (w->p->pt, &(w->p->maxfct), &(w->p->mnum), &(w->p->mtype),
             &(w->p->phase), &(w->p->n), &ddum, w->p->Ui, w->p->Uj, &idum, &nrhs,
             w->p->iparm, &(w->p->msglvl), &ddum, &ddum, 
             &(w->p->error), w->p->dparm);
  }

  if (w->p->Uj) PDOS_free(w->p->Uj);
  if (w->p->Ui) PDOS_free(w->p->Ui);
  if (w->p->U)  PDOS_free(w->p->U);
  if (w->p->tmp)PDOS_free(w->p->tmp);
  PDOS_free(w->p);
  w->p = NULL;
}

void formKKT(Work * w){
  /* ONLY UPPER TRIANGULAR PART IS STUFFED
   * forms column compressed KKT matrix
   * assumes column compressed form A matrix
   *
   * forms upper triangular part of [-I A'; A I]
   * 
   * puts into 1-based index'd *row compressed* form.
   */
  int j, k, kk;
  int nnzk;

  /* allocate memory */
  w->p->Unz = w->m + w->n + w->Ap[w->n];
  w->p->Ui = PDOS_malloc(sizeof(int)*w->p->n);
  w->p->Uj = PDOS_malloc(sizeof(int)*w->p->Unz);
  w->p->U  = PDOS_malloc(sizeof(double)*w->p->Unz);
  
  /* -I at top left and A^T at tope right: CCS to RCS */
  kk = 0;
  w->p->Ui[0] = 1;
  for(j = 0; j < w->n; j++) {
    nnzk = w->Ap[j+1] - w->Ap[j];
    w->p->Ui[j+1] = (1 + nnzk + w->p->Ui[j]);  /* nnz in row k */
    w->p->Uj[kk] = j+1;
    w->p->U[kk] = -1;
    kk++;

    /* copy entire column over */ 
    for(k = w->Ap[j]; k < w->Ap[j+1]; k++) {
      w->p->Uj[kk] = w->Ai[k] + w->n + 1; /* 1-based index */
      w->p->U[kk] = w->Ax[k];
      kk++;
    }
  }

  /* once here, we've copied first n rows */

  /* we copy the last m rows */
  /* I at bottom right */
  for(j = w->n; j < w->n + w->m; j++) {
    w->p->Ui[j+1] = (1 + w->p->Ui[j]);
    w->p->Uj[kk] = j+1;
    w->p->U[kk] = 1;
    kk++;
  }
 
#ifndef NDEBUG
  assert(w->p->Ui[w->p->n] - 1 == w->p->Unz);
  assert(kk == w->p->Unz);  
#endif
}


static inline void prepArgument(Work *w) {
  // uses the x and stilde memory space to store the argument, since we don't
  // need that memory space anymore

  // w->x = (-(x-lambda*c), b - (s+y))
  idxint i;
  for (i = 0; i < w->n; i++) {
    // set x = -(x - lambda*c)
    w->x[i] = w->lambda*w->c[i] - w->x[i];
  }
  for (i = 0; i < w->m; i++) {
    // set stilde = (b - (s + lambda*y))
    w->stilde[i] = w->b[i] - (w->s[i] + w->lambda*w->y[i]);
  }
}


void projectLinSys(Work * w){
  //static timer lin_timer;
  // this only modifies w->s = s + y - b
  //tic(&lin_timer);
  prepArgument(w);
  //toc(&lin_timer);

  //tic(&lin_timer);
  // solve the linear system
  static int idum;
  static const int nrhs = 1;

  pardiso (w->p->pt, &(w->p->maxfct), &(w->p->mnum), &(w->p->mtype), 
           &(w->p->phase), &(w->p->n), w->p->U, w->p->Ui, w->p->Uj,
           &idum, &nrhs, w->p->iparm, &(w->p->msglvl), w->x, w->p->tmp,
           &(w->p->error), w->p->dparm);

  if (w->p->error != 0 &&
      w->p->error != -100 &&  /* for hitting max iters */
      w->p->error != -101)
  {
    PDOS_printf("ERROR during solution: %d\n", w->p->error);
    exit(3);
  }
  //toc(&lin_timer);

  //tic(&lin_timer);
  /* now, stilde contains b - v - A*x,
   * where v = s + lambda * y. 
   * To recover stilde = b - A*x, we simply add "v"
   * to the solution.
   */
  // stilde = b - A*x
  static idxint i = 0;
  for (i = 0; i < w->m; i++) {
    // set stilde = (b - (s + lambda*y))
    w->stilde[i] += (w->s[i] + w->lambda*w->y[i]);
  }
  //toc(&lin_timer);
}

