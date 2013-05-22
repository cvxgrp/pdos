#include "mex.h" 
#include "matrix.h"
#include "pdos.h"

double getParameterField(const mxArray *params, const char *field, Data *d, Cone *k)
{
  const mxArray *tmp = mxGetField(params, 0, field);
  if(tmp == NULL) {
    mxFree(d); mxFree(k);
    mexErrMsgIdAndTxt("PDOS:getParams", "Params struct must contain a(n) `%s` entry.", field);
  }
  return *mxGetPr(tmp);
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
    mxFree(d); mxFree(k);
    mexErrMsgTxt("Data struct must contain a `A` entry.");
  }
  if (!mxIsSparse(A_mex)){
    mxFree(d); mxFree(k);
    mexErrMsgTxt("Input matrix A must be in sparse format (pass in sparse(A))");
  }
  const mxArray *b_mex = (mxArray *) mxGetField(data,0,"b");
  if(b_mex == NULL) {
    mxFree(d); mxFree(k);
    mexErrMsgTxt("Data struct must contain a `b` entry.");
  }
  const mxArray *c_mex = (mxArray *) mxGetField(data,0,"c"); 
  if(c_mex == NULL) {
    mxFree(d); mxFree(k);
    mexErrMsgTxt("Data struct must contain a `c` entry.");
  }

  
  const mxArray *cone = prhs[1];
  const mxArray *params = prhs[2];
  d->n = *(mxGetDimensions(c_mex));
  d->m = *(mxGetDimensions(b_mex));

  d->b = mxGetPr(b_mex);
  d->c = mxGetPr(c_mex);
  
  d->p->ALPHA = getParameterField(params, "ALPHA", d, k);
  d->p->MAX_ITERS = (idxint)getParameterField(params, "MAX_ITERS", d, k);
  
  d->p->EPS_ABS = getParameterField(params, "EPS_ABS", d, k);
  d->p->EPS_REL = getParameterField(params, "EPS_REL", d, k);

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
  k->qsize = *(mxGetDimensions(mxGetField(cone,0,"q")));
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
  
  plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
  mxSetPr(plhs[0], sol->x);
  mxSetM(plhs[0], d->n); 
  mxSetN(plhs[0], 1); 
  
  plhs[1] = mxCreateDoubleMatrix(0, 0, mxREAL);
  mxSetPr(plhs[1], sol->s);
  mxSetM(plhs[1], d->m); 
  mxSetN(plhs[1], 1); 

  plhs[2] = mxCreateDoubleMatrix(0, 0, mxREAL);
  mxSetPr(plhs[1], sol->y);
  mxSetM(plhs[1], d->m); 
  mxSetN(plhs[1], 1); 

  plhs[3] = mxCreateString(sol->status);
  
  mxFree(d->p); mxFree(d); mxFree(k->q); mxFree(k);
  
  //free(d->Ai);free(d->Ap);free(d);free(k->q);free(k);
  
  return; 
}

