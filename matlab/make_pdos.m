% the option -Dprintf=mexPrintf redfines all calls to "printf" with
%   "mexPrintf"
% switch(computer)
%     case {'PCWIN', 'GLNX86', 'MACI'}
%         evalc('system(''make -C .. purge packages CFLAGS="-m32 -DMATLAB_MEX_FILE"'')')
%     case {'PCWIN64', 'GLNXA64', 'SOL64', 'MACI64'}
%         evalc('system(''make -C .. purge packages CFLAGS="-DMATLAB_MEX_FILE"'')')
%     otherwise
%         evalc('system(''make -C .. purge packages CFLAGS="-DMATLAB_MEX_FILE"'')')
% end

%compile direct
d = '' ;
if (~isempty (strfind (computer, '64')))
    d = '-largeArrayDims -m64' ;
else
    d = '-m32';
end

i = sprintf ('-I../ -I../direct/external/SuiteSparse_config -I../direct/external/AMD/Include -I../direct/external/LDL/Include') ;
cmd = sprintf ('mex -v -O CFLAGS="-std=c99 -DDLONG -DMATLAB_MEX_FILE %s" %s', d, i) ;
amd_files = {'amd_order', 'amd_dump', 'amd_postorder', 'amd_post_tree', ...
    'amd_aat', 'amd_2', 'amd_1', 'amd_defaults', 'amd_control', ...
    'amd_info', 'amd_valid', 'amd_global', 'amd_preprocess' } ;
for i = 1 : length (amd_files)
    cmd = sprintf ('%s ../direct/external/AMD/Source/%s.c', cmd, amd_files {i}) ;
end
cmd = sprintf ('%s ../direct/external/LDL/Source/ldl.c  pdos_mex.c ../pdos.c ../cones.c ../cs.c ../linAlg.c ../util.c ../direct/private.c -lm -o pdos_direct', cmd) ;

eval (cmd) ;

% compile indirect
mex -v -O -largeArrayDims CFLAGS="\$CFLAGS -std=c99 -DMATLAB_MEX_FILE" pdos_mex.c ../pdos.c ../cones.c ../cs.c ../linAlg.c ../util.c ../indirect/private.c -I../ -o pdos_indirect -lm 
