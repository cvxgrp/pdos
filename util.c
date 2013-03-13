#include "util.h"

static clock_t tic_timestart;

void tic(void) {
  tic_timestart = clock();
}

float toc(void) {
  clock_t tic_timestop;
  tic_timestop = clock();
  printf("time: %8.4f seconds.\n", (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC);
  return (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC;
}

float tocq(void) {
  clock_t tic_timestop;
  tic_timestop = clock();
  return (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC;
}

void printConeData(Cone * k){
	int i;
	printf("num zeros = %i\n",k->f);
	printf("num LP = %i\n",k->l);
	printf("num SOCs = %i\n",k->qsize);
	printf("soc array:\n");
	for ( i=0;i<k->qsize;i++){
		printf("%i\n",k->q[i]);
	}
}

void printData(Data * d){
	printf("d->n is %i\n",d->n);
	printf("d->m is %i\n",d->m);
	printf("d->b[0] is %4f\n",d->b[0]);
	printf("d->c[0] is %4f\n",d->c[0]);
  printf("d->MAX_ITERS is %i\n",d->MAX_ITERS);
	printf("d->CG_MAX_ITS is %i\n",d->CG_MAX_ITS);
  printf("d->ALPH is %6f\n",d->ALPH);
  printf("d->EPS_ABS is %6f\n",d->EPS_ABS);
  //printf("d->EPS_ABS is %6f\n",d->EPS_ABS);
  //printf("d->UNDET_TOL is %6f\n",d->UNDET_TOL);
  printf("d->CG_TOL is %6f\n",d->CG_TOL);
  printf("d->Ap[0] is %i\n",d->Ap[0]);
  printf("d->Ap[1] is %i\n",d->Ap[1]);
	printf("d->Ai[0] is %i\n",d->Ai[0]);
	printf("d->Ax[0] is %4f\n",d->Ax[0]);
}

void printAll(Data * d, Work * w){
	int i;
	printf("\n z-half is \n");
	for( i=0;i<w->l;i++){
		printf("%f\n",w->z_half[i]);
	}  
	printf("\n z is \n");
	for( i=0;i<w->l;i++){
		printf("%f\n",w->z[i]);
	}
	printf("\n u is \n");
	for( i=0;i<w->l;i++){
		printf("%f\n",w->u[i]);
	}
}

