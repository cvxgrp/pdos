#include <Python.h>
#include "pdos.h"
#include "cvxopt.h"

// TODO: when normalizing, make a copy

/* WARNING: this code uses numpy array types
 *
 * WARNING: this code also does not check that the data for the matrix A is
 * actually column compressed storage for a sparse matrix. if it's not, the
 * code will just crash inelegantly. support for cvxopt matrix or scipy sparse
 * is planned, but not likely to be implemented soon.
 */


static inline void freeDataAndConeOnly(Data *d, Cone *k) {
  // this function is useful since the Data and Cone "structs" do not own the
  // memory for the arrays; numpy does.
  if(d) free(d);
  if(k->q) free(k->q);
  if(k) free(k);
  d = NULL; k = NULL;
}

// TODO: use PyObject * to keep track of whether or not two objects are equivalent (for warm-starting)
// static const double *prev_Ax, *prev_b, *prev_c;
// static const int *prev_Ai, *prev_Ap, *prev_q;

static Sol *solution = NULL;

static void cleanup()
{
  free_sol(solution);
}

static PyObject *solve(PyObject* self, PyObject *args, PyObject *keywords)
{
  /* Expects a function call 
   *     sol = solve(c,G,h,dims,opts)
   * The uppercase keywords are optional. If INDIRECT is #define'd, then
   * CG_MAX_ITS and CG_TOL can also be provided as options.
   *
   * `c` is a cvxopt (dense) column vector
   * `G` is a cvxopt (sparse) matrix
   * `h` is a cvxopt (dense) column vector
   * `dims` is a dictionary with
   *    `dims['f']` an integer giving the number of free variables
   *    `dims['l']` an integer specifying the dimension of positive orthant cone
   *    `dims['q']` an *array* specifying dimensions of second-order cones
   * `opts` is an optional dictionary with
   *    `opts['MAX_ITERS']` is an integer. Sets the maximum number of ADMM iterations.
   *        Defaults to 2000.
   *    `opts['EPS_ABS']` is a double. Sets the quitting tolerance for ADMM. 
   *        Defaults to 1e-3.
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
   *        Defaults to 1e-3.
   *  
   * The code returns a Python dictionary with three keys, 'x', 'y', and 'status'.
   * These report the primal and dual solution (as cvxopt dense matrices) and the solver
   * status (as a string).
   */
     
  Data *d = calloc(1,sizeof(Data)); // sets everything to 0
  Cone *k = calloc(1,sizeof(Cone)); // sets everything to 0
  // set default values
  d->MAX_ITERS = 2000;
  d->EPS_ABS = 1e-3;
  d->ALPH = 1.0;
  d->VERBOSE = 1;
#ifdef INDIRECT
  d->CG_MAX_ITS = 20;
  d->CG_TOL = 1e-3;
#endif
  
  matrix *c, *h;
  spmatrix *G;
  PyObject *dims, *opts = NULL;
  
  idxint m, n, i;
  
  if( !PyArg_ParseTuple(args, "OOOO!|O!",
        &c,
        &G,
        &h,
        &PyDict_Type, &dims,
        &PyDict_Type, &opts)
    ) { freeDataAndConeOnly(d,k); return NULL; }
        
  /* set G */
  if ((SpMatrix_Check(G) && SP_ID(G) != DOUBLE)){
      PyErr_SetString(PyExc_TypeError, "G must be a sparse 'd' matrix");
      freeDataAndConeOnly(d,k); return NULL;
  }
  if ((m = SP_NROWS(G)) <= 0) {
      PyErr_SetString(PyExc_ValueError, "m must be a positive integer");
      freeDataAndConeOnly(d,k); return NULL;
  }
  if ((n = SP_NCOLS(G)) <= 0) {
      PyErr_SetString(PyExc_ValueError, "n must be a positive integer");
      freeDataAndConeOnly(d,k); return NULL;
  }
  d->Ax = SP_VALD(G);
  d->Ai = SP_ROW(G);
  d->Ap = SP_COL(G);

  /* set c */
  if (!Matrix_Check(c) || MAT_NCOLS(c) != 1 || MAT_ID(c) != DOUBLE) {
      PyErr_SetString(PyExc_TypeError, "c must be a dense 'd' matrix with one column");
      freeDataAndConeOnly(d,k); return NULL;
  }

  if (MAT_NROWS(c) != n){
      PyErr_SetString(PyExc_ValueError, "c has incompatible dimension with G");
      freeDataAndConeOnly(d,k); return NULL;
  }
  d->c = MAT_BUFD(c);

  /* set h */
  if (!Matrix_Check(h) || MAT_NCOLS(h) != 1 || MAT_ID(h) != DOUBLE) {
    PyErr_SetString(PyExc_TypeError, "h must be a dense 'd' matrix with one column");
    freeDataAndConeOnly(d,k); return NULL;
  }

  if (MAT_NROWS(h) != m){
      PyErr_SetString(PyExc_ValueError, "h has incompatible dimension with G");
      freeDataAndConeOnly(d,k); return NULL;
  }
  d->b = MAT_BUFD(h);
  
  // set dimensions
  d->m = m; d->n = n;
  
  /* get dims['f'] */
  PyObject *freeObj = PyDict_GetItemString(dims, "f");
  if(freeObj) {
    if(PyInt_Check(freeObj) && ((k->f = (idxint) PyInt_AsLong(freeObj)) >= 0)) {
      // do nothing
    } else {
      PyErr_SetString(PyExc_TypeError, "dims['f'] ought to be a nonnegative integer");
      freeDataAndConeOnly(d,k); return NULL;
    }
  }
  
  /* get dims['l'] */
  PyObject *linearObj = PyDict_GetItemString(dims, "l");
  if(linearObj) {
    if(PyInt_Check(linearObj) && ((k->l = (idxint) PyInt_AsLong(linearObj)) >= 0)) {
      // do nothing
    } else {
      PyErr_SetString(PyExc_TypeError, "dims['l'] ought to be a nonnegative integer");
      freeDataAndConeOnly(d,k); return NULL;
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
            // do nothing
          } else {
            PyErr_SetString(PyExc_TypeError, "dims['q'] ought to be a list of positive integers");
            freeDataAndConeOnly(d,k); return NULL;
          }

      }
    } else {
      PyErr_SetString(PyExc_TypeError, "dims['q'] ought to be a list");
      freeDataAndConeOnly(d,k); return NULL;
    }
  }
  
  if(opts) {
    PyObject *dictObj = NULL;
    /* MAX_ITERS */
    dictObj = PyDict_GetItemString(dims, "MAX_ITERS");
    if(dictObj) {
      if(PyInt_Check(dictObj) && ((d->MAX_ITERS = (idxint) PyInt_AsLong(dictObj)) >= 0)) {
        // do nothing
      } else {
        PyErr_SetString(PyExc_TypeError, "opts['MAX_ITERS'] ought to be a nonnegative integer");
        freeDataAndConeOnly(d,k); return NULL;
      }
    }
    /* VERBOSE */
    dictObj = PyDict_GetItemString(dims, "VERBOSE");
    if(dictObj) {
      if(PyBool_Check(dictObj)) {
        d->VERBOSE = (idxint) PyInt_AsLong(dictObj);
      } else {
        PyErr_SetString(PyExc_TypeError, "opts['VERBOSE'] ought to be a boolean");
        freeDataAndConeOnly(d,k); return NULL;
      }
    }
    /* NORMALIZE */
    dictObj = PyDict_GetItemString(dims, "NORMALIZE");
    if(dictObj) {
      if(PyBool_Check(dictObj)) {
        d->NORMALIZE = (idxint) PyInt_AsLong(dictObj);
      } else {
        PyErr_SetString(PyExc_TypeError, "opts['NORMALIZE'] ought to be a boolean");
        freeDataAndConeOnly(d,k); return NULL;
      }
    }
    
    /* EPS_ABS */
    dictObj = PyDict_GetItemString(dims, "EPS_ABS");
    if(dictObj) {
      if(PyFloat_Check(dictObj) && ((d->EPS_ABS = (double) PyFloat_AsDouble(dictObj)) >= 0.0)) {
        // do nothing
      } else {
        PyErr_SetString(PyExc_TypeError, "opts['EPS_ABS'] ought to be a positive floating point value");
        freeDataAndConeOnly(d,k); return NULL;
      }
    }
    
    /* ALPHA */
    dictObj = PyDict_GetItemString(dims, "ALPHA");
    if(dictObj) {
      if(PyFloat_Check(dictObj)) {
        d->ALPH = (double) PyFloat_AsDouble(dictObj);
        if(d->ALPH >= 2.0 || d->ALPH <= 0.0) {
          PyErr_SetString(PyExc_TypeError, "opts['ALPHA'] ought to be a floating point value between 0 and 2 (noninclusive)");
          freeDataAndConeOnly(d,k); return NULL;
        }
      } else {
        PyErr_SetString(PyExc_TypeError, "opts['ALPHA'] ought to be a floating point value");
        freeDataAndConeOnly(d,k); return NULL;
      }
    }

#ifdef INDIRECT
    /* CG_MAX_ITS */
    dictObj = PyDict_GetItemString(dims, "CG_MAX_ITS");
    if(dictObj) {
      if(PyInt_Check(dictObj) && ((d->CG_MAX_ITS = (idxint) PyInt_AsLong(dictObj)) >= 0)) {
        // do nothing
      } else {
        PyErr_SetString(PyExc_TypeError, "opts['CG_MAX_ITS'] ought to be a nonnegative integer");
        freeDataAndConeOnly(d,k); return NULL;
      }
    }

    /* CG_TOL */
    dictObj = PyDict_GetItemString(dims, "CG_TOL");
    if(dictObj) {
      if(PyFloat_Check(dictObj) && ((d->CG_TOL = (double) PyFloat_AsDouble(dictObj)) >= 0.0)) {
        // do nothing
      } else {
        PyErr_SetString(PyExc_TypeError, "opts['CG_TOL'] ought to be a positive floating point value");
        freeDataAndConeOnly(d,k); return NULL;
      }
    }
#endif
  }
  
  // solve the problem
  // TODO: preserve the workspace
  solution = pdos(d, k);
  
  /* x */
  matrix *x;
  if(!(x = Matrix_New(n,1,DOUBLE)))
    return PyErr_NoMemory();
  for(i = 0; i < n; ++i)
    MAT_BUFD(x)[i] = solution->x[i];
        
  /* y */
  matrix *y;
  if(!(y = Matrix_New(m,1,DOUBLE)))
    return PyErr_NoMemory();
  for(i = 0; i < m; ++i)
    MAT_BUFD(y)[i] = solution->y[i];
  

  PyObject *returnDict = Py_BuildValue("{s:O,s:O,s:s}","x", x, "y", y, "status", solution->status);
  // give up ownership to the return dictionary
  Py_DECREF(x); Py_DECREF(y); 
  
  // do some cleanup
  freeDataAndConeOnly(d,k);

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
