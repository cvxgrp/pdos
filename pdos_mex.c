#include "mex.h" 
#include "matrix.h"
#include "pdos.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* matlab usage: ConeOS(A,b,c,cone,params); */
  /*
  if (nrhs != 4){
    mexErrMsgTxt("Four arguments are required in this order: A, b, c, cone struct");
  }
  */
  Data * d = malloc(sizeof(Data)); 
  Cone * k = malloc(sizeof(Cone));   
  const mxArray * A_mex = prhs[0];
  if (!mxIsSparse(A_mex)){
    mexErrMsgTxt("Input matrix A must be in sparse format (pass in sparse(A))");
  }
  const mxArray * b_mex = prhs[1];
  const mxArray * c_mex = prhs[2]; 
  const mxArray * cone = prhs[3];
  const mxArray * params = prhs[4];
  d->n = *(mxGetDimensions(c_mex));
  d->m = *(mxGetDimensions(b_mex));

  d->b = mxGetPr(b_mex);
  d->c = mxGetPr(c_mex);

  d->ALPH = (double)*mxGetPr(mxGetField(params,0,"ALPHA"));
  //d->UNDET_TOL = (double)*mxGetPr(mxGetField(params,0,"UNDET_TOL"));
  d->MAX_ITERS = (int)*mxGetPr(mxGetField(params,0,"MAX_ITERS"));
  d->EPS_ABS = (double)*mxGetPr(mxGetField(params,0,"EPS_ABS"));
  //d->EPS_REL = (double)*mxGetPr(mxGetField(params,0,"EPS_REL"));

  d->CG_MAX_ITS = (int)*mxGetPr(mxGetField(params,0,"CG_MAX_ITS"));
  d->CG_TOL = (double)*mxGetPr(mxGetField(params,0,"CG_TOL"));

  k->f = (int)*mxGetPr(mxGetField(cone,0,"f"));
  k->l = (int)*mxGetPr(mxGetField(cone,0,"l"));
  
  double * q_mex = mxGetPr(mxGetField(cone,0,"q"));
  k->qsize = *(mxGetDimensions(mxGetField(cone,0,"q")));
  int i;
  k->q = malloc(sizeof(int)*k->qsize);
  for ( i=0; i<k->qsize; i++ ){
    k->q[i] = (int)q_mex[i]; 
  }
  
  int Anz = mxGetNzmax(A_mex);
  d->Ax = (double *)mxGetPr(A_mex);
  d->Ap = (int *)malloc(sizeof(int)*Anz);
  d->Ai = (int *)malloc(sizeof(int)*Anz);
  long * A_i = (long*)mxGetIr(A_mex);
  /* XXX fix this crap: */
  for (i = 0; i < Anz; i++){
    d->Ai[i] = (int)A_i[i];
  }
  long * A_p = (long*)mxGetJc(A_mex);
  for (i = 0; i < d->n+1; i++) {          
    d->Ap[i] = (int)A_p[i];
  }

  /* printConeData(d,k); */
  /* printData(d); */
  Sol * sol = pdos(d,k);

  plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
  mxSetPr(plhs[0], sol->x);
  mxSetM(plhs[0], d->n); 
  mxSetN(plhs[0], 1); 

  plhs[1] = mxCreateDoubleMatrix(0, 0, mxREAL);
  mxSetPr(plhs[1], sol->z);
  mxSetM(plhs[1], d->m); 
  mxSetN(plhs[1], 1); 

  plhs[2] = mxCreateString(sol->status);
  
  free(d->Ai);free(d->Ap);free(d);free(k->q);free(k); 
  return; 
}

