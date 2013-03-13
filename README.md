Primal Dual Operator Splitting Method for Conic Optimization
============================================================
Brendan O'Donoghue, Eric Chu, and Stephen Boyd
----------------------------------------------

This code provides a solver for second-order cone problems. It is an
implementation of the algorithm described in [this
paper](http://www.stanford.edu/~boyd/). It provides both a direct and an
indirect solver in the form of a static library for inclusion in other
projects.

It simultaneously solves the primal cone program

    minimize     c'*x
    subject to   A*x + s == b
                 s in K 
                 
and its dual

    maximize     -b'*y
    subject to   -A'*y == c
                 y in K^* 

where `K` is a product cone of free cones, linear cones `{ x | x >= 0 }`, and second-order cones `{ (t,x) | ||x||_2 <= t }`; `K^*` is its dual cone.

Installing
----------
Typing `make` at the command line should do the trick. It will produce two libaries, `libpdosdir.a` and `libpdosindir.a`.

File an issue with us if the build process fails for you.

### Compiling a Matlab mex file
This is a little more difficult, depending on whether you have an older version of Matlab (32-bit only) or a newer version. For the newer versions, running `make_mex` (not yet included) will produce a usable mex file.

If `make_mex` fails and complains about an incompatible architecture... (pass in compiler flags...?)

Remember to include this directory in your Matlab path if you wish to use the mex file in your Matlab code.

### Compiling a Python extension
If you have no trouble with `make`, then this should be straightforward as well. The command

    python setup.py install

should install the two extensions, `pdos_direct` and `pdos_indirect` into your Python directories. You can run python and `import pdos_direct` or `import pdos_indirect` to use the codes. These two modules expose only a single `solve` function.

(XXX.... explain solve) 

Usage in C
----------
If `make` completes successfully, it will produce two static library files,
`libpdosdir.a` and `libpdosindir.a`. To include the libraries in your own
source code, compile with the linker option `-lpdosdir` or `-lpdosindir` (as
needed).

These libraries (and `pdos.h`) expose only three API functions:

* Sol *pdos(Data * d, Cone * k)
    
    This solves the problem specified in the `Data` and `Cone` structures. 
    The solution is returned in a `Sol` structure.
    
* void free_data(Data *d, Cone *k)
    
    This frees the `Data` and `Cone` structures.
    
* void free_sol(Sol *sol)

    This frees the `Sol` structure.
    
The three relevant data structures are:

    typedef struct PROBLEM_DATA {
      int n, m; /* problem dimensions */
      /* problem data, A, b, c: */
      double * Ax;
      int * Ai, * Ap;
      double * b, * c;
      /* problem parameters */
      int MAX_ITERS, CG_MAX_ITS;
      double EPS_ABS, EPS_INFEAS, ALPH, CG_TOL;
      int VERBOSE;
    } Data;

    typedef struct SOL_VARS {
      /* primal solution x, dual solution y */
      double * x, * y;
      char * status;
    } Sol;

    typedef struct Cone_t {
        int f;          /* number of linear equality constraints */
        int l;          /* length of LP cone */
        int *q;   		  /* array of second-order cone constraints */
        int qsize;      /* length of SOC array */
    } Cone;

The data matrix `A` is specified in column-compressed format and the vectors
`b` and `c` are specified as dense arrays. The solutions `x` (primal) and `y`
(dual) are returned as dense arrays. Cones are specified in terms of their
lengths; the only special one is the second-order cone, where the lengths are
provided as an array of second-order cone lengths (and a variable `qsize`
giving its length).


Scalability
-----------
Note that this code is merely meant as an implementation of the ideas in our
paper. The actual code does not use more than a single CPU. Nevertheless, for
problems that fit in memory on a single computer, this code will (attempt to)
solve them.

To scale this solver, one must either provide a distributed solver for linear
systems or a distributed matrix-vector multiplication.


TODO List
---------
Just a random list of items that we might do in the future....

* CVX shim?
* mex file interface?
* compile with -DINDIRECT to get an indirect solver
* multithreaded version?
* (Python) need some way of indicating warm start / saving state

Known Issues
------------
* When using OSX's built-in version of Python, `setup.py` may attempt to build Power PC versions of the Python extensions. If you upgraded to Xcode 4, then support for Power PC has been removed. In that case, the build process may complain about the Power PC architecture. Simply use the following instead:

    ARCHFLAGS='-arch i386 -arch x86_64' python setup.py install
