#include "util.h"

static clock_t tic_timestart;

void tic(void) {
  tic_timestart = clock();
}

float toc(void) {
  clock_t tic_timestop;
  tic_timestop = clock();
  PDOS_printf("time: %8.4f seconds.\n", (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC);
  return (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC;
}

float tocq(void) {
  clock_t tic_timestop;
  tic_timestop = clock();
  return (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC;
}

void printConeData(Cone * k){
	int i;
	PDOS_printf("num zeros = %i\n",k->f);
	PDOS_printf("num LP = %i\n",k->l);
	PDOS_printf("num SOCs = %i\n",k->qsize);
	PDOS_printf("soc array:\n");
	for ( i=0;i<k->qsize;i++){
		PDOS_printf("%i\n",k->q[i]);
	}
}

void printData(Data * d){
	PDOS_printf("d->n is %i\n",d->n);
	PDOS_printf("d->m is %i\n",d->m);
	PDOS_printf("d->b[0] is %4f\n",d->b[0]);
	PDOS_printf("d->c[0] is %4f\n",d->c[0]);
  PDOS_printf("d->MAX_ITERS is %i\n",d->MAX_ITERS);
	PDOS_printf("d->CG_MAX_ITS is %i\n",d->CG_MAX_ITS);
  PDOS_printf("d->ALPH is %6f\n",d->ALPH);
  PDOS_printf("d->EPS_ABS is %6f\n",d->EPS_ABS);
  //PDOS_printf("d->EPS_ABS is %6f\n",d->EPS_ABS);
  //PDOS_printf("d->UNDET_TOL is %6f\n",d->UNDET_TOL);
  PDOS_printf("d->CG_TOL is %6f\n",d->CG_TOL);
  PDOS_printf("d->Ap[0] is %i\n",d->Ap[0]);
  PDOS_printf("d->Ap[1] is %i\n",d->Ap[1]);
	PDOS_printf("d->Ai[0] is %i\n",d->Ai[0]);
	PDOS_printf("d->Ax[0] is %4f\n",d->Ax[0]);
}

void printAll(Data * d, Work * w){
	int i;
	PDOS_printf("\n z-half is \n");
	for( i=0;i<w->l;i++){
		PDOS_printf("%f\n",w->z_half[i]);
	}  
	PDOS_printf("\n z is \n");
	for( i=0;i<w->l;i++){
		PDOS_printf("%f\n",w->z[i]);
	}
	PDOS_printf("\n u is \n");
	for( i=0;i<w->l;i++){
		PDOS_printf("%f\n",w->u[i]);
	}
}

