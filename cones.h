#ifndef CONES_H_GUARD                                                              
#define CONES_H_GUARD

#include <string.h>
#include "globals.h"


typedef struct CONE {
    idxint f;          /* number of linear equality constraints */
    idxint l;          /* length of LP cone */
    idxint *q;   		   /* array of second-order cone constraints */
    idxint qsize;      /* length of SOC array */
} Cone;

void projCone(double *x,const Cone *k);

#endif
