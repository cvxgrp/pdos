function [data, cone, params] = read_pdos_data(name)

fi = fopen(name,'r');
n = fread(fi,1,                 'int64=>double');
m = fread(fi,1,                 'int64=>double');
cone.f = fread(fi,1,            'int64=>double');
cone.l = fread(fi,1,            'int64=>double');
k = fread(fi,1,                 'int64=>double');
cone.q = fread(fi,k,            'int64=>double');
data.b = fread(fi,m,            'double');
data.c = fread(fi,n,            'double');

params.MAX_ITERS = fread(fi,1,  'int64=>double');
params.CG_MAX_ITS = fread(fi,1, 'int64=>double');
params.ALPHA = fread(fi,1,      'double');
params.EPS_ABS = fread(fi,1,    'double');
params.CG_TOL = fread(fi,1,     'double');
params.VERBOSE = fread(fi,1,    'int64');
params.NORMALIZE = fread(fi,1,  'int64');

% triplet A
%{
[i,j,s] = find(data.A);
i = i-1;j=j-1;
NNZ=length(i);
fprintf(fi,'%u ',NNZ);
fprintf(fi,'\n');
fprintf(fi,'%u ',i');
fprintf(fi,'\n');
fprintf(fi,'%u ',j');
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',s');
fprintf(fi,'\n');
%}

% col-compressed A
NNZ = fread(fi, 1,  'int64=>double');
i = fread(fi, NNZ,  'int64=>double');    % length nnz
pw = fread(fi, n+1, 'int64=>double');    % length n+1
s = fread(fi, NNZ,  'double');  % lenght nnz

data.Ai = i;
data.Apw = pw;
data.As = s;


%{
[i,j,s] = find(W);
i = i-1;
tmp = full(sum(W)~=0);
pw = [0 cumsum(tmp)];
NNZ=length(i);
fprintf(fi,'%u ',NNZ);
fprintf(fi,'\n');
fprintf(fi,'%u ',i');
fprintf(fi,'\n');
fprintf(fi,'%u ',pw);
fprintf(fi,'\n');
fprintf(fi,'%6.18f ',s');
fprintf(fi,'\n');
%}
fclose(fi);
