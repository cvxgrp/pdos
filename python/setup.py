from distutils.core import setup, Extension
from glob import glob
from numpy import get_include

direct = Extension('pdos_direct', libraries = ['m'],
                    include_dirs = ['../', get_include()],
                    define_macros = [('PYTHON', None)],
                    sources = ['pdosmodule.c',
                        '../linAlg.c', '../cones.c', '../cs.c', 
                        '../common.c', '../pdos.c', '../util.c'
                    ] + glob('../direct/*.c'),
                    extra_compile_args=['-std=c99'])

indirect = Extension('pdos_indirect', libraries = ['m'],
                    include_dirs = ['../', get_include()],
                    define_macros = [('INDIRECT', None), ('PYTHON', None)],
                    sources = ['pdosmodule.c',
                        '../indirect/private.c',
                        '../linAlg.c', '../cones.c', '../cs.c', 
                        '../common.c', '../pdos.c', '../util.c'
                    ],
                    extra_compile_args=['-std=c99'])


setup(  name = 'pdos',
        version = '1.0',
        description = 'This is Python package to wrap our first-order solvers',
        ext_modules = [direct, indirect],
        requires = ["numpy (>= 1.7)"])
