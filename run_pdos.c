#include "pdos.h"
#include "run_pdos.h"

#ifndef DEMO_PATH
#define DEMO_PATH "../data_pdos"
#endif 

int main(int argc, char **argv)
{
  tic();
  FILE * fp;
  if(open_file(argc, argv, 1, DEMO_PATH, &fp)==-1) return -1;
  Cone * k = malloc(sizeof(Cone));
  Data * d = malloc(sizeof(Data));
  read_in_data(fp,d,k);
  fclose(fp);

  Sol * sol = pdos(d,k);
  printf("Total factorize + solve time %4f seconds\n",tocq());
  free_data(d,k);
  free_sol(sol);
  return 0;
}

void read_in_data(FILE * fp,Data * d, Cone * k){
#ifdef DLONG
	/* MATRIX IN DATA FILE MUST BE IN COLUMN COMPRESSED FORMAT */
	fscanf(fp, "%li", &(d->n));
	fscanf(fp, "%li", &(d->m));
  fscanf(fp, "%li", &(k->f));
  fscanf(fp, "%li", &(k->l));
  fscanf(fp, "%li", &(k->qsize)); 
  k->q = malloc(sizeof(idxint)*k->qsize);
  for(idxint i = 0; i < k->qsize; i++)
  { 
    fscanf(fp, "%li", &k->q[i]);
  }

  d->b = malloc(sizeof(double)*d->m);
  for(idxint i = 0; i < d->m; i++)
  { 
    fscanf(fp, "%lf", &d->b[i]);
  }

  d->c = malloc(sizeof(double)*d->n);
  for(idxint i = 0; i < d->n; i++)
  { 
    fscanf(fp, "%lf", &d->c[i]);
  }
  fscanf(fp, "%li", &(d->MAX_ITERS));
  fscanf(fp, "%li", &(d->CG_MAX_ITS)); 

  fscanf(fp, "%lf", &(d->ALPH));
  // fscanf(fp, "%lf", &(d->UNDET_TOL)); 
  fscanf(fp, "%lf", &(d->EPS_ABS)); 
  fscanf(fp, "%lf", &(d->EPS_INFEAS));
  fscanf(fp, "%lf", &(d->CG_TOL));
  fscanf(fp, "%li", &(d->VERBOSE));
  fscanf(fp, "%li", &(d->NORMALIZE));
    
  idxint Anz;
  fscanf(fp, "%li", &Anz);
	d->Ai = malloc(sizeof(idxint)*Anz);
  for(idxint i = 0; i < Anz; i++)
	{
		fscanf(fp, "%li", &d->Ai[i]);
	}
  d->Ap = malloc(sizeof(idxint)*(d->n+1));
	for(idxint i = 0; i < d->n+1; i++) 
	{
		fscanf(fp, "%li", &d->Ap[i]);
	}
  d->Ax = malloc(sizeof(double)*Anz);
	for(int i = 0; i < Anz; i++)
	{
		fscanf(fp, "%lf", &d->Ax[i]);
	}
#else
	/* MATRIX IN DATA FILE MUST BE IN COLUMN COMPRESSED FORMAT */
	fscanf(fp, "%i", &(d->n));
	fscanf(fp, "%i", &(d->m));
  fscanf(fp, "%i", &(k->f));
  fscanf(fp, "%i", &(k->l));
  fscanf(fp, "%i", &(k->qsize)); 
  k->q = malloc(sizeof(idxint)*k->qsize);
  for(int i = 0; i < k->qsize; i++)
  { 
    fscanf(fp, "%i", &k->q[i]);
  }

  d->b = malloc(sizeof(double)*d->m);
  for(int i = 0; i < d->m; i++)
  { 
    fscanf(fp, "%lf", &d->b[i]);
  }

  d->c = malloc(sizeof(double)*d->n);
  for(int i = 0; i < d->n; i++)
  { 
    fscanf(fp, "%lf", &d->c[i]);
  }
  fscanf(fp, "%i", &(d->MAX_ITERS));
  fscanf(fp, "%i", &(d->CG_MAX_ITS)); 

  fscanf(fp, "%lf", &(d->ALPH));
  // fscanf(fp, "%lf", &(d->UNDET_TOL)); 
  fscanf(fp, "%lf", &(d->EPS_ABS)); 
  fscanf(fp, "%lf", &(d->EPS_INFEAS));
  fscanf(fp, "%lf", &(d->CG_TOL));
  fscanf(fp, "%i", &(d->VERBOSE));
  fscanf(fp, "%i", &(d->NORMALIZE));
  
  int Anz;
  fscanf(fp, "%i", &Anz);
	d->Ai = malloc(sizeof(idxint)*Anz);
  for(int i = 0; i < Anz; i++)
	{
		fscanf(fp, "%i", &d->Ai[i]);
	}
  d->Ap = malloc(sizeof(idxint)*(d->n+1));
	for(int i = 0; i < d->n+1; i++) 
	{
		fscanf(fp, "%i", &d->Ap[i]);
	}
  d->Ax = malloc(sizeof(double)*Anz);
	for(int i = 0; i < Anz; i++)
	{
		fscanf(fp, "%lf", &d->Ax[i]);
	}
#endif
 
  /*	
  fscanf(fp, "%zu", &NNZ);
  int *Kr = malloc(sizeof(idxint)*NNZ);
	for(int i = 0; i < NNZ; i++)
	{
		fscanf(fp, "%i", &Kr[i]);
	}
	int *Kp=malloc(sizeof(idxint)*(w->l+1));
	for(int i = 0; i < w->l+1; i++)
	{
		fscanf(fp, "%i", &Kp[i]);
	}
	double *Kx=malloc(sizeof(double)*NNZ);
	for(int i = 0; i < NNZ; i++)
	{
		fscanf(fp, "%lf", &Kx[i]);
	}
  */
}

int open_file(int argc, char ** argv, int idx, char * default_file, FILE ** fb) 
{
	if (argc<idx+1){
		printf("Not enough arguments supplied, using %s as default\n", default_file);
	}
	else{
		*fb = fopen(argv[idx], "r");
		if (*fb != NULL) return 0;
		else{
			printf("Couldn't open file %s, using %s as default\n", argv[idx],default_file);
			fclose(*fb);
		}
	}
	*fb = fopen(default_file, "r");
	if (*fb == NULL){
		printf("Couldn't open %s\n",default_file);
		return -1;
	}
	return 0;
}
