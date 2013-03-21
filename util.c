#include "util.h"

static clock_t tic_timestart;

void tic(void) {
  tic_timestart = clock();
}

float toc(void) {
  clock_t tic_timestop;
  tic_timestop = clock();
  pdos_printf("time: %8.4f seconds.\n", (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC);
  return (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC;
}

float tocq(void) {
  clock_t tic_timestop;
  tic_timestop = clock();
  return (float)(tic_timestop - tic_timestart) / CLOCKS_PER_SEC;
}

void printConeData(Cone * k){
	int i;
	pdos_printf("num zeros = %i\n",k->f);
	pdos_printf("num LP = %i\n",k->l);
	pdos_printf("num SOCs = %i\n",k->qsize);
	pdos_printf("soc array:\n");
	for ( i=0;i<k->qsize;i++){
		pdos_printf("%i\n",k->q[i]);
	}
}

void printData(Data * d){
	pdos_printf("d->n is %i\n",d->n);
	pdos_printf("d->m is %i\n",d->m);
	pdos_printf("d->b[0] is %4f\n",d->b[0]);
	pdos_printf("d->c[0] is %4f\n",d->c[0]);
  pdos_printf("d->MAX_ITERS is %i\n",d->MAX_ITERS);
	pdos_printf("d->CG_MAX_ITS is %i\n",d->CG_MAX_ITS);
  pdos_printf("d->ALPH is %6f\n",d->ALPH);
  pdos_printf("d->EPS_ABS is %6f\n",d->EPS_ABS);
  //pdos_printf("d->EPS_ABS is %6f\n",d->EPS_ABS);
  //pdos_printf("d->UNDET_TOL is %6f\n",d->UNDET_TOL);
  pdos_printf("d->CG_TOL is %6f\n",d->CG_TOL);
  pdos_printf("d->Ap[0] is %i\n",d->Ap[0]);
  pdos_printf("d->Ap[1] is %i\n",d->Ap[1]);
	pdos_printf("d->Ai[0] is %i\n",d->Ai[0]);
	pdos_printf("d->Ax[0] is %4f\n",d->Ax[0]);
}

void printAll(Data * d, Work * w){
	int i;
	pdos_printf("\n z-half is \n");
	for( i=0;i<w->l;i++){
		pdos_printf("%f\n",w->z_half[i]);
	}  
	pdos_printf("\n z is \n");
	for( i=0;i<w->l;i++){
		pdos_printf("%f\n",w->z[i]);
	}
	pdos_printf("\n u is \n");
	for( i=0;i<w->l;i++){
		pdos_printf("%f\n",w->u[i]);
	}
}

