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
common_pdos = 'pdos_mex.c ../pdos.c ../cones.c ../cs.c ../util.c';


if (~isempty (strfind (computer, '64')))
    d = '-fPIC' ;
    arr = '-largeArrayDims';
else
    d = '-fPIC -m32';
    arr = '';
end

if ( isunix && ~ismac ) 
    link = '-lm -lrt';
else
    link = '-lm';
end

cmd = sprintf ('mex -v -O %s CFLAGS="-std=c99 -O3 -D_GNU_SOURCE -DMATLAB_MEX_FILE -DLDL_LONG -DDLONG %s" -I../', arr, d) ;
amd_files = {'amd_order', 'amd_dump', 'amd_postorder', 'amd_post_tree', ...
    'amd_aat', 'amd_2', 'amd_1', 'amd_defaults', 'amd_control', ...
    'amd_info', 'amd_valid', 'amd_global', 'amd_preprocess' } ;
for i = 1 : length (amd_files)
    cmd = sprintf ('%s ../direct/%s.c', cmd, amd_files {i}) ;
end

cmd = sprintf ('%s ../direct/ldl.c %s ../direct/private.c %s -o pdos_direct', cmd, common_pdos, link) ;

eval(cmd) ;

% compile indirect
cmd = sprintf('mex -v -O %s CFLAGS="-std=c99 -O3 -D_GNU_SOURCE -DMATLAB_MEX_FILE -DLDL_LONG -DDLONG %s" %s ../indirect/private.c -I../ -o pdos_indirect %s', arr, d, common_pdos, link);
eval(cmd);
