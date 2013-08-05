/* NB: this is a subset of the routines in the CSPARSE package by
  Tim Davis et. al., for the full package please visit
  http://www.cise.ufl.edu/research/sparse/CSparse/ */

/*
CSparse: a Concise Sparse matrix package.
Copyright (c) 2006, Timothy A. Davis.
http://www.suitesparse.com

--------------------------------------------------------------------------------

CSparse is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

CSparse is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this Module; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#ifndef CS_H_GUARD
#define CS_H_GUARD

#include "globals.h"

typedef struct cs_sparse    /* matrix in compressed-column or triplet form */
{
    idxint nzmax ;     /* maximum number of entries */
    idxint m ;         /* number of rows */
    idxint n ;         /* number of columns */
    idxint *p ;        /* column pointers (size n+1) or col indices (size nzmax) */
    idxint *i ;        /* row indices, size nzmax */
    double *x ;     /* numerical values, size nzmax */
    idxint nz ;        /* # of entries in triplet matrix, -1 for compressed-col */
} cs ;

cs *cs_compress (const cs *T);
cs *cs_done (cs *C, void *w, void *x, idxint ok);
cs *cs_spalloc (idxint m, idxint n, idxint nzmax, idxint values, idxint triplet);
cs *cs_spfree (cs *A);
double cs_cumsum (idxint *p, idxint *c, idxint n);
idxint *cs_pinv (idxint const *p, idxint n);
cs *cs_symperm (const cs *A, const idxint *pinv, idxint values);
cs *cs_transpose (const cs *A, idxint values);


#endif
