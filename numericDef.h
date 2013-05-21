#ifndef NUMERIC_H_GUARD                                                              
#define NUMERIC_H_GUARD

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
