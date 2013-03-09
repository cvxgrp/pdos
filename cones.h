#ifndef CONES_H_GUARD                                                              
#define CONES_H_GUARD

#include <string.h>

typedef struct Cone_t {
    int f;          /* number of linear equality constraints */
    int l;          /* length of LP cone */
    int *q;   		/* array of second-order cone constraints */
    int qsize;      /* length of SOC array */
} Cone;

void projCone(double * z,Cone* k);
void projDualCone(double * z,Cone* k);

#endif
