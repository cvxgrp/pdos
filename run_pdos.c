#include "pdos.h"
#include "run_pdos.h"

#ifndef DEMO_PATH
#define DEMO_PATH "../data/portfolio_test_01_10x100"
#endif
static timer PDOS_timer;

int main(int argc, char **argv)
{
  tic(&PDOS_timer);
  FILE * fp;
  if(open_file(argc, argv, 1, DEMO_PATH, &fp)==-1) return -1;
  Cone * k = malloc(sizeof(Cone));
  Data * d = malloc(sizeof(Data));
  read_in_data(fp,&d,&k);
  fclose(fp);
  printf("File IO %4f seconds\n",tocq(&PDOS_timer));

  tic(&PDOS_timer);
  Sol * sol = pdos(d,k);
  printf("Total factorize + solve time %4f seconds\n",tocq(&PDOS_timer));

  freeData(&d,&k);
  freeSol(&sol);
  return 0;
}

void read_in_data(FILE * fp, Data ** d, Cone ** k){
	/* MATRIX IN DATA FILE MUST BE IN COLUMN COMPRESSED FORMAT */
  (*d)->p = malloc(sizeof(Params));

  fread(&((*d)->n), sizeof(idxint), 1, fp);
	fread(&((*d)->m), sizeof(idxint), 1, fp);
	fread(&((*k)->f), sizeof(idxint), 1, fp);
	fread(&((*k)->l), sizeof(idxint), 1, fp);
	fread(&((*k)->qsize), sizeof(idxint), 1, fp);

  (*k)->q = malloc(sizeof(idxint)*(*k)->qsize);
  fread((*k)->q, sizeof(idxint), (*k)->qsize, fp);

  (*d)->b = malloc(sizeof(double)*(*d)->m);
  fread((*d)->b, sizeof(double), (*d)->m, fp);

  (*d)->c = malloc(sizeof(double)*(*d)->n);
  fread((*d)->c, sizeof(double), (*d)->n, fp);

  fread(&((*d)->p->MAX_ITERS), sizeof(idxint), 1, fp);
  fread(&((*d)->p->CG_MAX_ITS), sizeof(idxint), 1, fp);
  fread(&((*d)->p->ALPHA), sizeof(double), 1, fp);
  fread(&((*d)->p->EPS_ABS), sizeof(double), 1, fp);
  //fscanf(fp, READ_FLOAT, &((*d)->p->EPS_REL));
  fread(&((*d)->p->CG_TOL), sizeof(double), 1, fp);
  fread(&((*d)->p->VERBOSE), sizeof(idxint), 1, fp);
  fread(&((*d)->p->NORMALIZE), sizeof(idxint), 1, fp);

  idxint Anz;
  fread(&Anz, sizeof(idxint), 1, fp);

  (*d)->Ai = malloc(sizeof(idxint)*Anz);
  fread((*d)->Ai, sizeof(idxint), Anz, fp);

  (*d)->Ap = malloc(sizeof(idxint)*((*d)->n+1));
  fread((*d)->Ap, sizeof(idxint), (*d)->n+1, fp);

  (*d)->Ax = malloc(sizeof(double)*Anz);
  fread((*d)->Ax, sizeof(double), Anz, fp);
  
  (*d)->x = NULL; (*d)->y = NULL; (*d)->s = NULL;



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
