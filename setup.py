from distutils.core import setup, Extension
from glob import glob

direct = Extension('pdos_direct', libraries = ['m'],
                    include_dirs = ['.', 
                        'direct/external/AMD/Include', 
                        'direct/external/LDL/Include',
                        'direct/external/SuiteSparse_config'],
                    sources = ['pdosmodule.c',
                        'direct/private.c',
                        'direct/external/LDL/Source/ldl.c',
                        'cones.c', 'cs.c', 'pdos.c', 'util.c'
                    ] + glob('direct/external/AMD/Source/*.c'))

indirect = Extension('pdos_indirect', libraries = ['m'],
                    include_dirs = ['.'],
                    define_macros = [('INDIRECT', None)],
                    sources = ['pdosmodule.c',
                        'indirect/private.c',
                        'cones.c', 'cs.c', 'pdos.c', 'util.c'
                    ])


setup(  name = 'Primal-Dual Operator Splitting for Conic Programming',
        version = '1.0',
        description = 'This is Python package to wrap our first-order solvers',
        ext_modules = [direct, indirect])