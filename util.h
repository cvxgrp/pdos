#ifndef UTIL_H_GUARD
#define UTIL_H_GUARD

#include "globals.h"
#include <stdlib.h>
#include <stdio.h>

/* timing code courtesty of A. Domahidi */
#if (defined WIN32 || defined _WIN64 || defined _WINDLL)

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

/* Use POSIX clock_gettime() for timing on non-Windows machines */
#include <time.h>

typedef struct timer{
	struct timespec tic;
	struct timespec toc;
} timer;

#endif


void tic(timer *t); 
double toc(timer *t); 
double tocq(timer *t); 

#include "pdos.h"

void printConeData(const Cone * k);
void printData(const Data * d);
void printAll(const Data * d, const Work * w);


#endif
