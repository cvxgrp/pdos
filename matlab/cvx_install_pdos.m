cvx_clear
make_pdos
[ ver, isoctave, fs, ps ] = cvx_version;

copyfile('cvx_pdos.m',strcat(fs,'/shims'));
mkdir(strcat(fs,'/pdos'));
copyfile(strcat('pdos_direct.',ps),strcat(fs,'/pdos'));
copyfile(strcat('pdos_indirect.',ps),strcat(fs,'/pdos'));
cvx_setup

% mini test script
%{
cvx_begin
cvx_solver 'pdos'
variable x
minimize(x)
x>=1
cvx_end
%}
