from distutils.core import setup, Extension
from glob import glob

direct = Extension('pdos_direct', libraries = ['m'],
                    include_dirs = ['../'],
                    define_macros = [
                        ('DLONG', None), 
                        ('LDL_LONG', None),
                        ('PYTHON', None)],
                    sources = ['pdosmodule.c',
                        '../cones.c', '../cs.c', 
                        '../pdos.c', '../util.c'
                    ] + glob('../direct/*.c'),
                    extra_compile_args=['-std=c99'])

indirect = Extension('pdos_indirect', libraries = ['m'],
                    include_dirs = ['../'],
                    define_macros = [
                        ('DLONG', None), 
                        ('LDL_LONG', None), 
                        ('INDIRECT', None), 
                        ('PYTHON', None)],
                    sources = ['pdosmodule.c',
                        '../indirect/private.c',
                        '../cones.c', '../cs.c', 
                        '../pdos.c', '../util.c'
                    ],
                    extra_compile_args=['-std=c99'])


setup(  name = 'pdos',
        version = '1.0',
        description = 'This is Python package to wrap our first-order solvers',
        ext_modules = [direct, indirect],
        requires = ["cvxopt (>= 1.1.5)"])
