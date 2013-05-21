#include "pdos.h"
#include "run_pdos.h"

#ifndef DEMO_PATH
#define DEMO_PATH "../data_pdos"
#endif 

#ifdef DLONG
#define READ_INT "%li"
#else
#define READ_INT "%i"
#endif

#define READ_FLOAT "%lf"

int main(int argc, char **argv)
{
  tic();
  FILE * fp;
  if(open_file(argc, argv, 1, DEMO_PATH, &fp)==-1) return -1;
  Cone * k = malloc(sizeof(Cone));
  Data * d = malloc(sizeof(Data));
  read_in_data(fp,&d,&k);
  fclose(fp);

  Sol * sol = pdos(d,k);
  printf("Total factorize + solve time %4f seconds\n",tocq());
  freeData(&d,&k);
  freeSol(&sol);
  return 0;
}

void read_in_data(FILE * fp,Data ** d, Cone ** k){
	/* MATRIX IN DATA FILE MUST BE IN COLUMN COMPRESSED FORMAT */
  (*d)->p = malloc(sizeof(Params));
  
	fscanf(fp, READ_INT, &((*d)->n));
	fscanf(fp, READ_INT, &((*d)->m));
  fscanf(fp, READ_INT, &((*k)->f));
  fscanf(fp, READ_INT, &((*k)->l));
  fscanf(fp, READ_INT, &((*k)->qsize)); 
  (*k)->q = malloc(sizeof(idxint)*(*k)->qsize);
  for(idxint i = 0; i < (*k)->qsize; i++)
  { 
    fscanf(fp, READ_INT, &(*k)->q[i]);
  }

  (*d)->b = malloc(sizeof(double)*(*d)->m);
  for(idxint i = 0; i < (*d)->m; i++)
  { 
    fscanf(fp, READ_FLOAT, &(*d)->b[i]);
  }

  (*d)->c = malloc(sizeof(double)*(*d)->n);
  for(idxint i = 0; i < (*d)->n; i++)
  { 
    fscanf(fp, READ_FLOAT, &(*d)->c[i]);
  }
  fscanf(fp, READ_INT, &((*d)->p->MAX_ITERS));
  fscanf(fp, READ_INT, &((*d)->p->CG_MAX_ITS)); 
  fscanf(fp, READ_INT, &((*d)->p->SEARCH_ITERS)); 

  fscanf(fp, READ_FLOAT, &((*d)->p->ALPHA));
  fscanf(fp, READ_FLOAT, &((*d)->p->BETA));
  fscanf(fp, READ_FLOAT, &((*d)->p->TAU));
  fscanf(fp, READ_FLOAT, &((*d)->p->EPS_ABS)); 
  //fscanf(fp, READ_FLOAT, &(d->EPS_INFEAS));
  fscanf(fp, READ_FLOAT, &((*d)->p->CG_TOL));
  fscanf(fp, READ_INT, &((*d)->p->VERBOSE));
  fscanf(fp, READ_INT, &((*d)->p->NORMALIZE));
    
  idxint Anz;
  fscanf(fp, READ_INT, &Anz);
	(*d)->Ai = malloc(sizeof(idxint)*Anz);
  for(idxint i = 0; i < Anz; i++)
	{
		fscanf(fp, READ_INT, &(*d)->Ai[i]);
	}
  (*d)->Ap = malloc(sizeof(idxint)*((*d)->n+1));
	for(idxint i = 0; i < (*d)->n+1; i++) 
	{
		fscanf(fp, READ_INT, &(*d)->Ap[i]);
	}
  (*d)->Ax = malloc(sizeof(double)*Anz);
	for(int i = 0; i < Anz; i++)
	{
		fscanf(fp, READ_FLOAT, &(*d)->Ax[i]);
	}

 
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
