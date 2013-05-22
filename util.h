#ifndef UTIL_H_GUARD
#define UTIL_H_GUARD

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "pdos.h"

void tic(void); 
double toc(void); 
double tocq(void); 
void printConeData(const Cone * k);
void printData(const Data * d);
void printAll(const Data * d, const Work * w);


#endif
