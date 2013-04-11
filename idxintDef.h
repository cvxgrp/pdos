#ifndef IDXINT_H_GUARD                                                              
#define IDXINT_H_GUARD

//#include "pdos.h"
#ifdef DLONG

#ifdef _WIN64
#define idxint __int64
#else
#define idxint long
#endif

#else
#define idxint int
#endif

#endif
