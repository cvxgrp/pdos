switch(computer)
    case {'PCWIN', 'GLNX86', 'MACI'}
        evalc('system(''make -C .. purge packages CFLAGS=-m32'')')
    case {'PCWIN64', 'GLNXA64', 'SOL64', 'MACI64'}
        evalc('system(''make -C .. purge packages'')')
    otherwise
        evalc('system(''make -C .. purge packages'')')
end

mex -v -O -largeArrayDims pdos_mex.c ../pdos.c ../cones.c ../cs.c ../linAlg.c ../util.c ../direct/private.c -I../ -I../direct/external/SuiteSparse_config -I../direct/external/AMD/Include -I../direct/external/LDL/Include -o pdos_direct -lm ../direct/external/AMD/Lib/libamd.a ../direct/external/LDL/Lib/libldl.a 

mex -v -O -largeArrayDims pdos_mex.c ../pdos.c ../cones.c ../cs.c ../linAlg.c ../util.c ../indirect/private.c -I../ -o pdos_indirect -lm 
