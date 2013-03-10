#ifndef UTIL_H_GUARD
#define UTIL_H_GUARD

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "pdos.h"

void tic(void); 
float toc(void); 
float tocq(void); 
void printConeData(Data * d,Cone * k);
void printData(Data * d);
void printAll(Data * d, Work * w);


#endif
