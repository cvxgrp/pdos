#include <Python.h>
#include "pdos.h"

/* WARNING: this code converts Python integers (which are stored as C long's)
 * into C ints. we don't anticipate working with arrays that push int
 * precision, but if that happens someday, we'll have to fix this.
 *
 * WARNING: this code also does not check that the data for the matrix A is
 * actually column compressed storage for a sparse matrix. if it's not, the
 * code will just crash inelegantly. support for cvxopt matrix or scipy sparse
 * is planned, but not likely to be implemented soon.
 */

static PyObject *PDOSError;

// converts a list of python floats into a double *
static int convertDoubleArray(PyObject *obj, void *data)
{
  int i, len;
  PyObject *tmp;
  double **p = (double **) data;
  if( !PySequence_Check(obj) ) {
    PyErr_SetString(PDOSError, "argument must be a sequence"); 
    return 0; 
  }
  
  len = (int) PySequence_Size(obj);
  if( len < 0 ) {
    PyErr_SetString(PDOSError, "unknown list length"); 
    return 0; 
  }
  
  // allocate memory
  *p = malloc(len*sizeof(double));
  
  // deep copy data
  for(i = 0; i < len; ++i) {
    tmp = PySequence_GetItem(obj, i); // creates a new reference
    if( !tmp ) {
      PyErr_SetString(PDOSError, "unable to access list");
      // free(*p);  // freed externally (when we call free_data)
      return 0; 
    }
    
    if( !PyFloat_Check(tmp) ) {
      PyErr_SetString(PDOSError, "must be a list of floats");
      // free(*p);  // freed externally (when we call free_data)
      return 0; 
    }
    
    (*p)[i] = PyFloat_AS_DOUBLE(tmp);
    Py_DECREF(tmp);  // decrement the reference
  }
  return 1;
}

// converts a list of python integers into an int *
static int convertIntegerArray(PyObject *obj, void *data)
{
  int i, len;
  PyObject *tmp;
  int **p = (int **) data;
  if( !PySequence_Check(obj) ) {
    PyErr_SetString(PDOSError, "argument must be a sequence"); 
    return 0; 
  }
  
  len = (int) PySequence_Size(obj);
  if( len < 0 ) {
    PyErr_SetString(PDOSError, "unknown list length"); 
    return 0; 
  }
  
  // allocate memory
  *p = malloc(len*sizeof(int));
  
  // deep copy data
  for(i = 0; i < len; ++i) {
    tmp = PySequence_GetItem(obj, i); // creates a new reference
    if( !tmp ) {
      PyErr_SetString(PDOSError, "unable to access list");
      // free(*p);  // freed externally (when we call free_data)
      return 0; 
    }
    
    if( !PyInt_Check(tmp) ) {
      PyErr_SetString(PDOSError, "must be a list of ints");
      // free(*p);  // freed externally (when we call free_data)
      return 0; 
    }
    
    // WARNING: python stores integers as longs, this eliminates precision
    // which is fine for arrays with integer values < 2^31.
    (*p)[i] = (int) PyInt_AS_LONG(tmp);
    Py_DECREF(tmp);  // decrement the reference
  }
  return 1;
}

// converts a list of python floats into a double * and stores "data.n"
static int convertPrimalObj(PyObject *obj, void *data)
{
  Data *d = (Data *) data;
  int return_val = convertDoubleArray(obj, &(d->c));
  if(return_val)
    d->n = (int) PySequence_Size(obj); // this has already succeeded
  
  return return_val;
}

// converts a list of python floats into a double * and stores "data.m"
static int convertDualObj(PyObject *obj, void *data)
{
  Data *d = (Data *) data;
  int return_val = convertDoubleArray(obj, &(d->b));
  if(return_val)
    d->m = (int) PySequence_Size(obj); // this has already succeeded

  return return_val;
}

// converts a list of python ints into an int* and stores "k.qsize"
static int convertCone(PyObject *obj, void *cone)
{
  Cone *k = (Cone *) cone;
  int return_val = convertIntegerArray(obj, &(k->q));
  if(return_val)
    k->qsize = (int) PySequence_Size(obj); // this has already succeeded

  return return_val;
}


static PyObject *solve(PyObject* self, PyObject *args, PyObject *keywords)
{
  /* Expects a function call 
   *     solve(Ax, Ai, Ap, b, c, f=,l=,q=, MAX_ITERS=, EPS_ABS=, EPS_INFEAS=, ALPHA=)
   * The uppercase keywords are optional. If INDIRECT is #define'd, then
   * CG_MAX_ITS and CG_TOL are also optional keyword arguments.
   *
   * "A" is a sparse matrix in column compressed storage. "Ax" are the values,
   * "Ai" are the rows, and "Ap" are the column pointers.
   * Ax is a list of doubles
   * Ai is a list of ints
   * Ap is a list of ints
   *
   * b is a (dense) list of doubles
   * c is a (dense) list of doubles
   *
   * f is an integer giving the number of free variables
   * l is an integer giving the number of nonnegative constraints
   * q is a list of integers giving the number of cone constraints
   * 
   * MAX_ITERS is an integer. Sets the maximum number of ADMM iterations.
   *  Defaults to 2000.
   * EPS_ABS is a double. Sets the quitting tolerance for ADMM. 
   *  Defaults to 1e-4.
   * EPS_INFEAS is a double. Sets the quitting tolerance for infeasibility.
   *  Defaults to 5e-5.
   * ALPHA is a double in (0,2) (non-inclusive). Sets the over-relaxation
   *  parameter. Defaults to 1.0.
   *
   * CG_MAX_ITS is an integer. Sets the maximum number of CG iterations.
   *  Defaults to 20.
   * CG_TOL is a double. Sets the tolerance for CG.
   *  Defaults to 1e-3.
   *  
   */
     
     
#ifdef INDIRECT
  static char *kwlist[] = {"Ax","Ai","Ap","b","c","f","l","q","MAX_ITERS", "EPS_ABS", "EPS_INFEAS", "ALPHA", "CG_MAX_ITS", "CG_TOL", NULL};
#else
  static char *kwlist[] = {"Ax","Ai","Ap","b","c","f","l","q","MAX_ITERS", "EPS_ABS", "EPS_INFEAS", "ALPHA", NULL};
#endif
  Data *d = calloc(1,sizeof(Data)); // sets everything to 0
  Cone *k = calloc(1,sizeof(Cone)); // sets everything to 0
  d->MAX_ITERS = 2000;
  d->EPS_ABS = 1e-4;
  d->EPS_INFEAS = 5e-5;
  d->ALPH = 1.0;
#ifdef INDIRECT
  d->CG_MAX_ITS = 20;
  d->CG_TOL = 1e-3;
#endif

#ifdef INDIRECT
  if( !PyArg_ParseTupleAndKeywords(args, keywords, "O&O&O&O&O&|iiO&idddid", kwlist,
      &convertDoubleArray, &(d->Ax),
      &convertIntegerArray, &(d->Ai),
      &convertIntegerArray, &(d->Ap),
      &convertPrimalObj, d,
      &convertDualObj, d,
      &(k->f),
      &(k->l),
      &convertCone, k,
      &(d->MAX_ITERS), 
      &(d->EPS_ABS),
      &(d->EPS_INFEAS),
      &(d->ALPH),
      &(d->CG_MAX_ITS),
      &(d->CG_TOL))
    ) { free_data(d,k); return NULL; }
  
#else
  if( !PyArg_ParseTupleAndKeywords(args, keywords, "O&O&O&O&O&|iiO&iddd", kwlist,
      &convertDoubleArray, &(d->Ax),
      &convertIntegerArray, &(d->Ai),
      &convertIntegerArray, &(d->Ap),
      &convertPrimalObj, d,
      &convertDualObj, d,
      &(k->f),
      &(k->l),
      &convertCone, k,
      &(d->MAX_ITERS), 
      &(d->EPS_ABS),
      &(d->EPS_INFEAS),
      &(d->ALPH))
    ) { free_data(d,k); return NULL; }
#endif
  // TODO: check that parameter values are correct
  Sol *tmp = pdos(d, k);
  
  printData(d);
  printConeData(k);

  if(tmp) {
    // TODO: build a new return object in python
    free_sol(tmp);
  }
  free_data(d,k);
  
  Py_INCREF(Py_None);

  return Py_None;
}

static PyMethodDef PDOSMethods[] =
{
  {"solve", (PyCFunction)solve, METH_VARARGS | METH_KEYWORDS, 
    "Solve a conic optimization problem."},
  {NULL, NULL, 0, NULL} // sentinel
};

PyMODINIT_FUNC
#ifdef INDIRECT
  initpdos_indirect(void)
#else
  initpdos_direct(void)
#endif  
{
  PyObject *m;

#ifdef INDIRECT
  m = Py_InitModule("pdos_indirect", PDOSMethods);
#else
  m = Py_InitModule("pdos_direct", PDOSMethods);
#endif
  
  if(m == NULL)
    return;

#ifdef INDIRECT
  PDOSError = PyErr_NewException("pdos_indirect.error", NULL, NULL);
#else
  PDOSError = PyErr_NewException("pdos_direct.error", NULL, NULL);
#endif

  Py_INCREF(PDOSError);
  PyModule_AddObject(m, "error", PDOSError);

}
