#include "util.h"

// timing code courtesty of A. Domahidi
#if (defined WIN32 || defined _WIN64)

/* Use Windows QueryPerformanceCounter for timing */
#include <Windows.h>

typedef struct timer{
	LARGE_INTEGER tic;
	LARGE_INTEGER toc;
	LARGE_INTEGER freq;
} timer;


#elif (defined __APPLE__)

#include <mach/mach_time.h>

/* Use MAC OSX  mach_time for timing */
typedef struct timer{
	uint64_t tic;
	uint64_t toc;
	mach_timebase_info_data_t tinfo;
} timer;

#else

/* Use POSIX clocl_gettime() for timing on non-Windows machines */
#include <time.h>

typedef struct timer{
	struct timespec tic;
	struct timespec toc;
} timer;

#endif

static timer PDOS_timer;

#if (defined WIN32 || _WIN64)

void tic()
{
	QueryPerformanceFrequency(&PDOS_timer.freq);
	QueryPerformanceCounter(&PDOS_timer.tic);
}

double tocq()
{
	QueryPerformanceCounter(&PDOS_timer.toc);
	return ((PDOS_timer.toc.QuadPart - PDOS_timer.tic.QuadPart) / (double)PDOS_timer.freq.QuadPart);
}

#elif (defined __APPLE__)

void tic()
{
    /* read current clock cycles */
    PDOS_timer.tic = mach_absolute_time();
}

double tocq()
{

    uint64_t duration; /* elapsed time in clock cycles*/
    
    PDOS_timer.toc = mach_absolute_time();
    duration = PDOS_timer.toc - PDOS_timer.tic;
    
    /*conversion from clock cycles to nanoseconds*/
    mach_timebase_info(&(PDOS_timer.tinfo));
    duration *= PDOS_timer.tinfo.numer;
    duration /= PDOS_timer.tinfo.denom;

    return (double)duration / 1000000000;
}


#else

/* read current time */
void tic()
{
	clock_gettime(CLOCK_MONOTONIC, &PDOS_timer.tic);
}


/* return time passed since last call to tic on this timer */
double tocq()
{
	struct timespec temp;
    
	clock_gettime(CLOCK_MONOTONIC, &PDOS_timer.toc);	
    
	if ((PDOS_timer.toc.tv_nsec - PDOS_timer.tic.tv_nsec)<0) {
		temp.tv_sec = PDOS_timer.toc.tv_sec - PDOS_timer.tic.tv_sec-1;
		temp.tv_nsec = 1000000000+PDOS_timer.toc.tv_nsec - PDOS_timer.tic.tv_nsec;
	} else {
		temp.tv_sec = PDOS_timer.toc.tv_sec - PDOS_timer.tic.tv_sec;
		temp.tv_nsec = PDOS_timer.toc.tv_nsec - PDOS_timer.tic.tv_nsec;
	}
	return (double)temp.tv_sec + (double)temp.tv_nsec / 1000000000;
}

#endif

double toc()
{
  double val = tocq();
  PDOS_printf("time: %8.4f seconds.\n", val);
  return val;
}

void printConeData(const Cone * k){
	idxint i;
#ifdef DLONG
	PDOS_printf("num zeros = %li\n",k->f);
	PDOS_printf("num LP = %li\n",k->l);
	PDOS_printf("num SOCs = %li\n",k->qsize);
	PDOS_printf("soc array:\n");
	for ( i=0;i<k->qsize;i++){
		PDOS_printf("%li\n",k->q[i]);
	}
#else
	PDOS_printf("num zeros = %i\n",k->f);
	PDOS_printf("num LP = %i\n",k->l);
	PDOS_printf("num SOCs = %i\n",k->qsize);
	PDOS_printf("soc array:\n");
	for ( i=0;i<k->qsize;i++){
		PDOS_printf("%i\n",k->q[i]);
	}
#endif
}

void printData(const Data * d){
#ifdef DLONG
	PDOS_printf("d->n is %li\n",d->n);
	PDOS_printf("d->m is %li\n",d->m);
	PDOS_printf("d->b[0] is %4f\n",d->b[0]);
	PDOS_printf("d->c[0] is %4f\n",d->c[0]);
  PDOS_printf("d->MAX_ITERS is %li\n",d->p->MAX_ITERS);
	PDOS_printf("d->CG_MAX_ITS is %li\n",d->p->CG_MAX_ITS);
  PDOS_printf("d->ALPH is %6f\n",d->p->ALPHA);
  PDOS_printf("d->EPS_ABS is %6f\n",d->p->EPS_ABS);
  //PDOS_printf("d->EPS_ABS is %6f\n",d->EPS_ABS);
  //PDOS_printf("d->UNDET_TOL is %6f\n",d->UNDET_TOL);
  PDOS_printf("d->CG_TOL is %6f\n",d->p->CG_TOL);
  PDOS_printf("d->Ap[0] is %li\n",d->Ap[0]);
  PDOS_printf("d->Ap[1] is %li\n",d->Ap[1]);
	PDOS_printf("d->Ai[0] is %li\n",d->Ai[0]);
	PDOS_printf("d->Ax[0] is %4f\n",d->Ax[0]);
#else
	PDOS_printf("d->n is %i\n",d->n);
	PDOS_printf("d->m is %i\n",d->m);
	PDOS_printf("d->b[0] is %4f\n",d->b[0]);
	PDOS_printf("d->c[0] is %4f\n",d->c[0]);
  PDOS_printf("d->MAX_ITERS is %i\n",d->p->MAX_ITERS);
	PDOS_printf("d->CG_MAX_ITS is %i\n",d->p->CG_MAX_ITS);
  PDOS_printf("d->ALPH is %6f\n",d->p->ALPHA);
  PDOS_printf("d->EPS_ABS is %6f\n",d->p->EPS_ABS);
  //PDOS_printf("d->EPS_ABS is %6f\n",d->EPS_ABS);
  //PDOS_printf("d->UNDET_TOL is %6f\n",d->UNDET_TOL);
  PDOS_printf("d->CG_TOL is %6f\n",d->p->CG_TOL);
  PDOS_printf("d->Ap[0] is %i\n",d->Ap[0]);
  PDOS_printf("d->Ap[1] is %i\n",d->Ap[1]);
	PDOS_printf("d->Ai[0] is %i\n",d->Ai[0]);
	PDOS_printf("d->Ax[0] is %4f\n",d->Ax[0]);
#endif
}

void printAll(const Data * d, const Work * w){
	idxint i;
	PDOS_printf("\n x is \n");
	for( i=0;i<d->n;i++){
		PDOS_printf("%f\n",w->x[i]);
	}  
	PDOS_printf("\n s is \n");
	for( i=0;i<d->m;i++){
		PDOS_printf("%f\n",w->s[i]);
	}
	PDOS_printf("\n y is \n");
	for( i=0;i<d->m;i++){
		PDOS_printf("%f\n",w->y[i]);
	}
}

