#ifndef GLOBALS_H_GUARD                                                              
#define GLOBALS_H_GUARD

// uncomment the following line to print KKT factorization information
// #define PRINTKKT

// redefine printfs and memory allocators as needed
#ifdef MATLAB_MEX_FILE
  #include "mex.h"
  #define PDOS_printf   mexPrintf
  #define PDOS_free     mxFree
  #define PDOS_malloc   mxMalloc
  #define PDOS_calloc   mxCalloc
#elif defined PYTHON
  #include <Python.h>
  #include <stdlib.h>
  #define PDOS_printf   PySys_WriteStdout
  #define PDOS_free     free
  #define PDOS_malloc   malloc
  #define PDOS_calloc   calloc
#else
  #include <stdio.h>
  #include <stdlib.h>
  #define PDOS_printf   printf
  #define PDOS_free     free
  #define PDOS_malloc   malloc
  #define PDOS_calloc   calloc
#endif

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
