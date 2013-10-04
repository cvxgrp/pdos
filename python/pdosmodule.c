#include <Python.h>
#include "pdos.h"
#include "cvxopt.h"

/* this code uses cvxopt matrix types */


static inline void freeDataAndConeOnly(Data **d, Cone **k) {
  // this function is useful since the Data and Cone "structs" do not own the
  // memory for the arrays; numpy does.
  if((*d)->p) free((*d)->p); 
  if(*d) free(*d);
  if((*k)->q) free((*k)->q);
  if(*k) free(*k);
  *d = NULL; *k = NULL;
}

// TODO: use PyObject * to keep track of whether or not two objects are equivalent (for warm-starting)

static Sol *solution = NULL;

static void cleanup()
{
  freeSol(&solution);
}

static PyObject *solve(PyObject* self, PyObject *args, PyObject *keywords)
{
  /* Expects a function call 
   *     sol = solve(c,G,h,dims,opts)
   * The uppercase keywords are optional. If INDIRECT is #define'd, then
   * CG_MAX_ITS and CG_TOL can also be provided as options.
   *
   * `c` is a cvxopt (dense) column vector
   * `A` is a cvxopt (sparse) matrix
   * `b` is a cvxopt (dense) column vector
   * `dims` is a dictionary with
   *    `dims['f']` an integer giving the number of equality constraints
   *    `dims['l']` an integer specifying the dimension of positive orthant cone
   *    `dims['q']` an *array* specifying dimensions of second-order cones
   * `opts` is an optional dictionary with
   *    `opts['MAX_ITERS']` is an integer. Sets the maximum number of ADMM iterations.
   *        Defaults to 2000.
   *    `opts['EPS_ABS']` is a double. Sets the quitting tolerance for ADMM. 
   *        Defaults to 5e-3.
   *    `opts['ALPHA']` is a double in (0,2) (non-inclusive). Sets the over-relaxation
   *        parameter. Defaults to 1.0.
   *    `opts['VERBOSE']` is an integer (or Boolean) either 0 or 1. Sets the verbosity of
   *        the solver. Defaults to 1 (or True).
   *    `opts['NORMALIZE']` is an integer (or Boolean) either 0 or 1. Tells the solver to
   *        normalize the data. Defaults to 0 (or False).
   *
   *    `opts['CG_MAX_ITS']` is an integer. Sets the maximum number of CG iterations.
   *        Defaults to 20.
   *    `opts['CG_TOL']` is a double. Sets the tolerance for CG.
   *        Defaults to 1e-4.
   *  
   * The code returns a Python dictionary with three keys, 'x', 'y', and 'status'.
   * These report the primal and dual solution (as cvxopt dense matrices) and the solver
   * status (as a string).
   */
     
  Data *d = calloc(1,sizeof(Data)); // sets everything to 0
  Cone *k = calloc(1,sizeof(Cone)); // sets everything to 0
  d->p = malloc(sizeof(Params));
  // set default values
  d->p->MAX_ITERS = 2000;
  d->p->EPS_ABS = 5e-3;
  d->p->ALPHA = 1.0;
  d->p->VERBOSE = 1;
#ifdef INDIRECT
  d->p->CG_MAX_ITS = 20;
  d->p->CG_TOL = 1e-4;
#endif
  
  matrix *c, *b;
  matrix *x0 = NULL, *y0 = NULL, *s0 = NULL;
  spmatrix *A;
  PyObject *dims, *opts = NULL;
  
  idxint m, n, i, num_conic_variables = 0;
  
  if( !PyArg_ParseTuple(args, "OOOO!|O!OOO",
        &c,
        &A,
        &b,
        &PyDict_Type, &dims,
        &PyDict_Type, &opts,
        &x0,
        &y0,
        &s0)
    ) { freeDataAndConeOnly(&d,&k); return NULL; }
        
  /* set A */
  if ((SpMatrix_Check(A) && SP_ID(A) != DOUBLE)){
      PyErr_SetString(PyExc_TypeError, "A must be a sparse 'd' matrix");
      freeDataAndConeOnly(&d,&k); return NULL;
  }
  if ((m = SP_NROWS(A)) <= 0) {
      PyErr_SetString(PyExc_ValueError, "m must be a positive integer");
      freeDataAndConeOnly(&d,&k); return NULL;
  }
  if ((n = SP_NCOLS(A)) <= 0) {
      PyErr_SetString(PyExc_ValueError, "n must be a positive integer");
      freeDataAndConeOnly(&d,&k); return NULL;
  }
  d->Ax = SP_VALD(A);
  d->Ai = SP_ROW(A);
  d->Ap = SP_COL(A);

  /* set c */
  if (!Matrix_Check(c) || MAT_NCOLS(c) != 1 || MAT_ID(c) != DOUBLE) {
      PyErr_SetString(PyExc_TypeError, "c must be a dense 'd' matrix with one column");
      freeDataAndConeOnly(&d,&k); return NULL;
  }

  if (MAT_NROWS(c) != n){
      PyErr_SetString(PyExc_ValueError, "c has incompatible dimension with A");
      freeDataAndConeOnly(&d,&k); return NULL;
  }
  d->c = MAT_BUFD(c);

  /* set h */
  if (!Matrix_Check(b) || MAT_NCOLS(b) != 1 || MAT_ID(b) != DOUBLE) {
    PyErr_SetString(PyExc_TypeError, "b must be a dense 'd' matrix with one column");
    freeDataAndConeOnly(&d,&k); return NULL;
  }

  if (MAT_NROWS(b) != m){
      PyErr_SetString(PyExc_ValueError, "b has incompatible dimension with A");
      freeDataAndConeOnly(&d,&k); return NULL;
  }
  d->b = MAT_BUFD(b);
  
  // set dimensions
  d->m = m; d->n = n;
  
  /* get dims['f'] */
  PyObject *freeObj = PyDict_GetItemString(dims, "f");
  if(freeObj) {
    if(PyInt_Check(freeObj) && ((k->f = (idxint) PyInt_AsLong(freeObj)) >= 0)) {
      num_conic_variables += k->f;
    } else {
      PyErr_SetString(PyExc_TypeError, "dims['f'] ought to be a nonnegative integer");
      freeDataAndConeOnly(&d,&k); return NULL;
    }
  }
  
  /* get dims['l'] */
  PyObject *linearObj = PyDict_GetItemString(dims, "l");
  if(linearObj) {
    if(PyInt_Check(linearObj) && ((k->l = (idxint) PyInt_AsLong(linearObj)) >= 0)) {
      num_conic_variables += k->l;
    } else {
      PyErr_SetString(PyExc_TypeError, "dims['l'] ought to be a nonnegative integer");
      freeDataAndConeOnly(&d,&k); return NULL;
    }
  }
  
  /* get dims['q'] */
  PyObject *socObj = PyDict_GetItemString(dims, "q");
  if(socObj) {
    if (PyList_Check(socObj)) {
      k->qsize = PyList_Size(socObj);
      k->q = calloc(k->qsize, sizeof(idxint));
      for (i = 0; i < k->qsize; ++i) {
          PyObject *qi = PyList_GetItem(socObj, i);
          if(PyInt_Check(qi) && ((k->q[i] = (idxint) PyInt_AsLong(qi)) > 0)) {
            num_conic_variables += k->q[i];
          } else {
            PyErr_SetString(PyExc_TypeError, "dims['q'] ought to be a list of positive integers");
            freeDataAndConeOnly(&d,&k); return NULL;
          }

      }
    } else {
      PyErr_SetString(PyExc_TypeError, "dims['q'] ought to be a list");
      freeDataAndConeOnly(&d,&k); return NULL;
    }
  }
  
  if( num_conic_variables != m ){
      PyErr_SetString(PyExc_ValueError, "Number of rows of G does not match dims.f+dims.l+sum(dims.q)");
      freeDataAndConeOnly(&d,&k); return NULL;
  }
  
  if(opts) {
    PyObject *dictObj = NULL;
    /* MAX_ITERS */
    dictObj = PyDict_GetItemString(opts, "MAX_ITERS");
    if(dictObj) {
      if(PyInt_Check(dictObj) && ((d->p->MAX_ITERS = (idxint) PyInt_AsLong(dictObj)) >= 0)) {
        // do nothing
      } else {
        PyErr_SetString(PyExc_TypeError, "opts['MAX_ITERS'] ought to be a nonnegative integer");
        freeDataAndConeOnly(&d,&k); return NULL;
      }
    }
    
    /* VERBOSE */
    dictObj = PyDict_GetItemString(opts, "VERBOSE");
    if(dictObj) {
      if(PyBool_Check(dictObj)) {
        d->p->VERBOSE = (idxint) PyInt_AsLong(dictObj);
      } else {
        PyErr_SetString(PyExc_TypeError, "opts['VERBOSE'] ought to be a boolean");
        freeDataAndConeOnly(&d,&k); return NULL;
      }
    }
    /* NORMALIZE */
    dictObj = PyDict_GetItemString(opts, "NORMALIZE");
    if(dictObj) {
      if(PyBool_Check(dictObj)) {
        d->p->NORMALIZE = (idxint) PyInt_AsLong(dictObj);
      } else {
        PyErr_SetString(PyExc_TypeError, "opts['NORMALIZE'] ought to be a boolean");
        freeDataAndConeOnly(&d,&k); return NULL;
      }
    }
    
    /* EPS_ABS */
    dictObj = PyDict_GetItemString(opts, "EPS_ABS");
    if(dictObj) {
      if(PyFloat_Check(dictObj) && ((d->p->EPS_ABS = (double) PyFloat_AsDouble(dictObj)) >= 0.0)) {
        // do nothing
      } else {
        PyErr_SetString(PyExc_TypeError, "opts['EPS_ABS'] ought to be a positive floating point value");
        freeDataAndConeOnly(&d,&k); return NULL;
      }
    }
    
    /* ALPHA */
    dictObj = PyDict_GetItemString(opts, "ALPHA");
    if(dictObj) {
      if(PyFloat_Check(dictObj)) {
        d->p->ALPHA = (double) PyFloat_AsDouble(dictObj);
        if(d->p->ALPHA >= 2.0 || d->p->ALPHA <= 0.0) {
          PyErr_SetString(PyExc_TypeError, "opts['ALPHA'] ought to be a floating point value between 0 and 2 (noninclusive)");
          freeDataAndConeOnly(&d,&k); return NULL;
        }
      } else {
        PyErr_SetString(PyExc_TypeError, "opts['ALPHA'] ought to be a floating point value");
        freeDataAndConeOnly(&d,&k); return NULL;
      }
    }

#ifdef INDIRECT
    /* CG_MAX_ITS */
    dictObj = PyDict_GetItemString(opts, "CG_MAX_ITS");
    if(dictObj) {
      if(PyInt_Check(dictObj) && ((d->p->CG_MAX_ITS = (idxint) PyInt_AsLong(dictObj)) >= 0)) {
        // do nothing
      } else {
        PyErr_SetString(PyExc_TypeError, "opts['CG_MAX_ITS'] ought to be a nonnegative integer");
        freeDataAndConeOnly(&d,&k); return NULL;
      }
    }

    /* CG_TOL */
    dictObj = PyDict_GetItemString(opts, "CG_TOL");
    if(dictObj) {
      if(PyFloat_Check(dictObj) && ((d->p->CG_TOL = (double) PyFloat_AsDouble(dictObj)) >= 0.0)) {
        // do nothing
      } else {
        PyErr_SetString(PyExc_TypeError, "opts['CG_TOL'] ought to be a positive floating point value");
        freeDataAndConeOnly(&d,&k); return NULL;
      }
    }
#endif
  }
  
  d->x = NULL; d->y = NULL; d->s = NULL;
  if(x0) { 
    /* set x0 */
    if (!Matrix_Check(x0) || MAT_NCOLS(x0) != 1 || MAT_ID(x0) != DOUBLE) {
        PyErr_SetString(PyExc_TypeError, "x0 must be a dense 'd' matrix with one column");
        freeDataAndConeOnly(&d,&k); return NULL;
    }

    if (MAT_NROWS(x0) != n){
        PyErr_SetString(PyExc_ValueError, "x0 has incompatible dimension with A");
        freeDataAndConeOnly(&d,&k); return NULL;
    }
    d->x = MAT_BUFD(x0);
  }
  if(y0) { 
    /* set y0 */
    if (!Matrix_Check(y0) || MAT_NCOLS(y0) != 1 || MAT_ID(y0) != DOUBLE) {
        PyErr_SetString(PyExc_TypeError, "y0 must be a dense 'd' matrix with one column");
        freeDataAndConeOnly(&d,&k); return NULL;
    }

    if (MAT_NROWS(y0) != m){
        PyErr_SetString(PyExc_ValueError, "y0 has incompatible dimension with A");
        freeDataAndConeOnly(&d,&k); return NULL;
    }
    d->y = MAT_BUFD(y0);
  }
  if(s0) { 
    /* set s0 */
    if (!Matrix_Check(s0) || MAT_NCOLS(s0) != 1 || MAT_ID(s0) != DOUBLE) {
        PyErr_SetString(PyExc_TypeError, "s0 must be a dense 'd' matrix with one column");
        freeDataAndConeOnly(&d,&k); return NULL;
    }

    if (MAT_NROWS(s0) != m){
        PyErr_SetString(PyExc_ValueError, "s0 has incompatible dimension with A");
        freeDataAndConeOnly(&d,&k); return NULL;
    }
    d->s = MAT_BUFD(s0);
  }
  
  // solve the problem
  // TODO: preserve the workspace
  solution = pdos(d, k);
  
  /* x */
  matrix *x;
  if(!(x = Matrix_New(n,1,DOUBLE)))
    return PyErr_NoMemory();
  memcpy(MAT_BUFD(x), solution->x, n*sizeof(double));
  
  /* s */
  matrix *s;
  if(!(s = Matrix_New(m,1,DOUBLE)))
    return PyErr_NoMemory();
  memcpy(MAT_BUFD(s), solution->s, m*sizeof(double));
  
        
  /* y */
  matrix *y;
  if(!(y = Matrix_New(m,1,DOUBLE)))
    return PyErr_NoMemory();
  memcpy(MAT_BUFD(y), solution->y, m*sizeof(double));

  PyObject *returnDict = Py_BuildValue("{s:O,s:O,s:O,s:s}","x", x, "s", s, "y", y, "status", solution->status);
  // give up ownership to the return dictionary
  Py_DECREF(x); Py_DECREF(s); Py_DECREF(y); 
  
  // do some cleanup
  freeDataAndConeOnly(&d,&k);

  return returnDict;
}

static PyMethodDef PDOSMethods[] =
{
  {"solve", (PyCFunction)solve, METH_VARARGS, 
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
  
  if (import_cvxopt() < 0) return; // for cvxopt support

  if(m == NULL)
    return;
  
  Py_AtExit(&cleanup);
}
