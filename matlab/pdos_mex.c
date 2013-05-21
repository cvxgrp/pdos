#include "mex.h" 
#include "matrix.h"
#include "pdos.h"

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
    free(d); free(k);
    mexErrMsgTxt("Data struct must contain a `c` entry.");
  }

  const mxArray *cone = prhs[1];
  const mxArray *params = prhs[2];
  d->n = *(mxGetDimensions(c_mex));
  d->m = *(mxGetDimensions(b_mex));

  d->b = mxGetPr(b_mex);
  d->c = mxGetPr(c_mex);
  
  d->p->ALPHA = (double)*mxGetPr(mxGetField(params,0,"ALPHA"));
  d->p->BETA = (double)*mxGetPr(mxGetField(params,0,"BETA"));
  d->p->TAU = (double)*mxGetPr(mxGetField(params,0,"TAU"));
  d->p->SEARCH_ITERS = (idxint)*mxGetPr(mxGetField(params,0,"SEARCH_ITERS"));
  //d->UNDET_TOL = (double)*mxGetPr(mxGetField(params,0,"UNDET_TOL"));
  d->p->MAX_ITERS = (idxint)*mxGetPr(mxGetField(params,0,"MAX_ITERS"));
  d->p->EPS_ABS = (double)*mxGetPr(mxGetField(params,0,"EPS_ABS"));
  //d->EPS_INFEAS = (double)*mxGetPr(mxGetField(params,0,"EPS_INFEAS"));

  d->p->CG_MAX_ITS = (idxint)*mxGetPr(mxGetField(params,0,"CG_MAX_ITS"));
  d->p->CG_TOL = (double)*mxGetPr(mxGetField(params,0,"CG_TOL"));
  d->p->VERBOSE = (idxint)*mxGetPr(mxGetField(params,0,"VERBOSE"));
  d->p->NORMALIZE = (idxint)*mxGetPr(mxGetField(params,0,"NORMALIZE"));

  k->f = (idxint)*mxGetPr(mxGetField(cone,0,"f"));
  k->l = (idxint)*mxGetPr(mxGetField(cone,0,"l"));
  
  double * q_mex = mxGetPr(mxGetField(cone,0,"q"));
  k->qsize = *(mxGetDimensions(mxGetField(cone,0,"q")));
  idxint i;
  k->q = mxMalloc(sizeof(idxint)*k->qsize);
  for ( i=0; i < k->qsize; i++ ){
    k->q[i] = (idxint)q_mex[i]; 
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
  mxSetPr(plhs[1], sol->y);
  mxSetM(plhs[1], d->m); 
  mxSetN(plhs[1], 1); 

  plhs[2] = mxCreateString(sol->status);
  
  mxFree(d->p); mxFree(d); mxFree(k->q); mxFree(k);
  
  //free(d->Ai);free(d->Ap);free(d);free(k->q);free(k);
  return; 
}

