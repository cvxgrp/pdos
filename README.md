Primal Dual Operator Splitting Method for Conic Optimization
============================================================

Installing
==========
* compile with -DINDIRECT to get an indirect solver
* Compiles into 64 bit.
* for 32 bit matlab, need to be a little more "careful"

Using as Library
================
* exposes "pdos" in a .a, static library file



* running make should do everything and produce a direct and an indirect solver
* run a demo to make sure it works
* no demo source code at the moment

* tested on MacOSX 10.8 (Mountain Lion)
* tested on Debian Linux

TODO:
* CVX shim?
* indirect solver
* matlab source code
* python interface