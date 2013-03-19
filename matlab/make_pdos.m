% the option -Dprintf=mexPrintf redfines all calls to "printf" with
%   "mexPrintf"
switch(computer)
    case {'PCWIN', 'GLNX86', 'MACI'}
        evalc('system(''make -C .. purge packages CFLAGS="-m32 -Dprintf=mexPrintf"'')')
    case {'PCWIN64', 'GLNXA64', 'SOL64', 'MACI64'}
        evalc('system(''make -C .. purge packages CFLAGS="-Dprintf=mexPrintf"'')')
    otherwise
        evalc('system(''make -C .. purge packages CFLAGS="-Dprintf=mexPrintf"'')')
end

% the option -std=c99 is to ensure we use the right c version
% the option -Dprintf=mexPrintf redfines all calls to "printf" with
%   "mexPrintf"
mex -v -O -largeArrayDims CFLAGS="\$CFLAGS -std=c99 -Dprintf=mexPrintf" pdos_mex.c ../pdos.c ../cones.c ../cs.c ../linAlg.c ../util.c ../direct/private.c -I../ -I../direct/external/SuiteSparse_config -I../direct/external/AMD/Include -I../direct/external/LDL/Include -o pdos_direct -lm ../direct/external/AMD/Lib/libamd.a ../direct/external/LDL/Lib/libldl.a 

mex -v -O -largeArrayDims CFLAGS="\$CFLAGS -std=c99 -Dprintf=mexPrintf" pdos_mex.c ../pdos.c ../cones.c ../cs.c ../linAlg.c ../util.c ../indirect/private.c -I../ -o pdos_indirect -lm 
