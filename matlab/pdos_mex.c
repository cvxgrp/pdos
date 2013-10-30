#include "mex.h" 
#include "matrix.h"
#include "pdos.h"

static double getParameterField(const mxArray *params, const char *field, Data *d, Cone *k)
{
  // helper function for getting a field of the params struct
  const mxArray *tmp = mxGetField(params, 0, field);
  if(tmp == NULL) {
    mxFree(d); mxFree(k);
    mexErrMsgIdAndTxt("PDOS:getParams", "Params struct must contain a(n) `%s` entry.", field);
  }
  return *mxGetPr(tmp);
}

static idxint getVectorLength(const mxArray *vec, const char *vec_name) {
  // get the vector length (even if transposed)
  long int m = mxGetM(vec);
  long int n = mxGetN(vec);
  
  if(m == 0 || n == 0) {
    return 0;
  }
  if(m == 1) {
    return n;
  }
  if (n == 1) {
    return m;
  }
  
  mexErrMsgIdAndTxt("PDOS:getVector", "Expected row or column vector `%s`.", vec_name);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* matlab usage: pdos(data,cone,params); */
  
  if (nrhs != 3){
    mexErrMsgTxt("Three arguments are required in this order: data struct, cone struct, params struct");
  }
  
  Data * d = mxMalloc(sizeof(Data)); 
  Cone * k = mxMalloc(sizeof(Cone));
  d->p = mxMalloc(sizeof(Params));
  const mxArray *data = prhs[0];
   
  const mxArray *A_mex = (mxArray *) mxGetField(data,0,"A");
  if(A_mex == NULL) {
    mxFree(d->p); mxFree(d); mxFree(k);
    mexErrMsgTxt("Data struct must contain a `A` entry.");
  }
  if (!mxIsSparse(A_mex)){
    mxFree(d->p); mxFree(d); mxFree(k);
    mexErrMsgTxt("Input matrix A must be in sparse format (pass in sparse(A))");
  }
  if (mxIsComplex(A_mex)) {
    mxFree(d->p); mxFree(d); mxFree(k);
    mexErrMsgTxt("Input matrix A cannot be complex");
  }
  
  const mxArray *b_mex = (mxArray *) mxGetField(data,0,"b");
  if(b_mex == NULL) {
    mxFree(d->p); mxFree(d); mxFree(k);
    mexErrMsgTxt("Data struct must contain a `b` entry.");
  }
  if(mxIsSparse(b_mex)) {
    mxFree(d->p); mxFree(d); mxFree(k);
    mexErrMsgTxt("Input vector b must be dense (pass in full(b))");
  }
  if (mxIsComplex(b_mex)) {
    mxFree(d->p); mxFree(d); mxFree(k);
    mexErrMsgTxt("Input vector b cannot be complex");
  }
  
  const mxArray *c_mex = (mxArray *) mxGetField(data,0,"c"); 
  if(c_mex == NULL) {
    mxFree(d->p); mxFree(d); mxFree(k);
    mexErrMsgTxt("Data struct must contain a `c` entry.");
  }
  if(mxIsSparse(c_mex)) {
    mxFree(d->p); mxFree(d); mxFree(k);
    mexErrMsgTxt("Input vector c must be dense (pass in full(c))");
  }
  if (mxIsComplex(c_mex)) {
    mxFree(d->p); mxFree(d); mxFree(k);
    mexErrMsgTxt("Input vector c cannot be complex");
  }
  
  const mxArray *x0 = (mxArray *) mxGetField(data,0,"x0");   
  if(x0 != NULL && mxIsSparse(x0)) {
    mxFree(d->p); mxFree(d); mxFree(k);
    mexErrMsgTxt("Initial vector x0 must be dense (pass in full(x0))");
  }
  if (x0 != NULL && mxIsComplex(x0)) {
    mxFree(d->p); mxFree(d); mxFree(k);
    mexErrMsgTxt("Initial vector x0 cannot be complex");
  }
  
  const mxArray *y0 = (mxArray *) mxGetField(data,0,"y0"); 
  if(y0 != NULL && mxIsSparse(y0)) {
    mxFree(d->p); mxFree(d); mxFree(k);
    mexErrMsgTxt("Initial vector y0 must be dense (pass in full(y0))");
  }
  if (y0 != NULL && mxIsComplex(y0)) {
    mxFree(d); mxFree(k);
    mexErrMsgTxt("Initial vector y0 cannot be complex");
  }

  
  const mxArray *s0 = (mxArray *) mxGetField(data,0,"s0"); 
  if(s0 != NULL && mxIsSparse(s0)) {
    mxFree(d); mxFree(k);
    mexErrMsgTxt("Initial vector s0 must be dense (pass in full(s0))");
  }
  if (s0 != NULL && mxIsComplex(s0)) {
    mxFree(d); mxFree(k);
    mexErrMsgTxt("Initial vector s0 cannot be complex");
  }
  
  const mxArray *cone = prhs[1];
  const mxArray *params = prhs[2];
  
  d->n = getVectorLength(c_mex,"data.c");
  d->m = getVectorLength(b_mex,"data.b");

  d->b = mxGetPr(b_mex);
  d->c = mxGetPr(c_mex);
  
  d->x = x0 ? mxGetPr(x0) : NULL;
  d->y = y0 ? mxGetPr(y0) : NULL;
  d->s = s0 ? mxGetPr(s0) : NULL;
  
  d->p->ALPHA = getParameterField(params, "ALPHA", d, k);
  d->p->MAX_ITERS = (idxint)getParameterField(params, "MAX_ITERS", d, k);
  
  d->p->EPS_ABS = getParameterField(params, "EPS_ABS", d, k);

  d->p->CG_MAX_ITS = (idxint)getParameterField(params, "CG_MAX_ITS", d, k);
  d->p->CG_TOL = getParameterField(params, "CG_TOL", d, k);
  d->p->VERBOSE = (idxint)getParameterField(params, "VERBOSE", d, k);
  d->p->NORMALIZE = (idxint)getParameterField(params, "NORMALIZE", d, k);


  const mxArray *f_mex = (mxArray *) mxGetField(cone,0,"f"); 
  if(f_mex == NULL) {
    mxFree(d); mxFree(k);
    mexErrMsgTxt("Cone struct must contain a `f` entry.");
  }
  k->f = (idxint)*mxGetPr(f_mex);
  const mxArray *l_mex = (mxArray *) mxGetField(cone,0,"l"); 
  if(l_mex == NULL) {
    mxFree(d); mxFree(k);
    mexErrMsgTxt("Cone struct must contain a `l` entry.");
  }
  k->l = (idxint)*mxGetPr(l_mex);
  
  const mxArray *q_mex = (mxArray *) mxGetField(cone,0,"q"); 
  if(q_mex == NULL) {
    mxFree(d); mxFree(k);
    mexErrMsgTxt("Cone struct must contain a `q` entry.");
  }
  
  double * q_mex_vals = mxGetPr(q_mex);
  k->qsize = getVectorLength(mxGetField(cone,0,"q"), "cone.q");
  idxint i;
  
  k->q = mxMalloc(sizeof(idxint)*k->qsize);
  for ( i=0; i < k->qsize; i++ ){
    k->q[i] = (idxint)q_mex_vals[i]; 
  }
  
  int Anz = mxGetNzmax(A_mex);
  d->Ax = (double *)mxGetPr(A_mex);
  d->Ai = (long*)mxGetIr(A_mex);
  d->Ap = (long*)mxGetJc(A_mex);

  /* printConeData(d,k); */
  /* printData(d); */
  
  Sol *sol = pdos(d,k);
  
  plhs[0] = mxCreateDoubleMatrix(d->n, 1, mxREAL);
  mxSetPr(plhs[0], sol->x);
  
  plhs[1] = mxCreateDoubleMatrix(d->m, 1, mxREAL);
  mxSetPr(plhs[1], sol->s);

  plhs[2] = mxCreateDoubleMatrix(d->m, 1, mxREAL);
  mxSetPr(plhs[2], sol->y);

  plhs[3] = mxCreateString(sol->status);
  
  mxFree(d->p); mxFree(d); mxFree(k->q); mxFree(k);
    
  return; 
}

