function write_pdos_data(data,cone,params,name)

MAX_ITERS = 2000; % maximum num iterations for admm
EPS_ABS   = 5e-4; % quitting tolerances
alpha=1.8;        % relaxation parameter (alpha = 1 is unrelaxed)
NORMALIZE = 1;

% conjugate gradient (CG) settings:
USE_INDIRECT = false; % use conjugate gradient rather than direct method
CG_MAX_ITS = 20; % max iterations for CG
CG_TOL = 1e-4; % max CG quitting tolerance
CG_VERBOSE = false; % CG prints summary
%%
if ~isfield(params,'MAX_ITERS');params.MAX_ITERS = MAX_ITERS;end
if ~isfield(params,'EPS_ABS');params.EPS_ABS = EPS_ABS;end
if ~isfield(params,'ALPHA');params.ALPHA = alpha;end
if ~isfield(params,'NORMALIZE');params.NORMALIZE = NORMALIZE;end
% CG:
if ~isfield(params,'USE_INDIRECT');params.USE_INDIRECT = USE_INDIRECT;end
if ~isfield(params,'CG_MAX_ITS');params.CG_MAX_ITS = CG_MAX_ITS;end
if ~isfield(params,'CG_TOL');params.CG_TOL = CG_TOL;end
if ~isfield(params,'CG_VERBOSE');params.CG_VERBOSE = CG_VERBOSE;end

n = length(data.c);
m = size(data.A,1);

%Q=sparse([zeros(n) data.A' data.c;
%    -data.A zeros(m,m) data.b;
%    -data.c' -data.b' 0]);
%W=sparse([speye(n+m+1) Q';Q -speye(n+m+1)]);

delete(name);
fi = fopen(name,'w');
fwrite(fi,n,                'int64');
fwrite(fi,m,                'int64');
fwrite(fi,cone.f,           'int64');
fwrite(fi,cone.l,           'int64');
fwrite(fi,length(cone.q),   'int64');
fwrite(fi,cone.q,           'int64');
fwrite(fi,data.b,           'double');
fwrite(fi,data.c,           'double');
fwrite(fi,params.MAX_ITERS, 'int64');
fwrite(fi,params.CG_MAX_ITS,'int64');
fwrite(fi,params.ALPHA,     'double');
fwrite(fi,params.EPS_ABS,   'double');
fwrite(fi,params.CG_TOL,    'double');
fwrite(fi,params.VERBOSE,   'int64');
fwrite(fi,params.NORMALIZE, 'int64');

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
[i,j,s] = find(sparse(data.A));
i = i-1;
tmp = full(sum(data.A~=0));
pw = [0 cumsum(tmp)];
NNZ=length(i);
fwrite(fi, NNZ, 'int64');
fwrite(fi, i,   'int64');    % length nnz
fwrite(fi, pw,  'int64');    % length n+1
fwrite(fi, s,   'double');  % lenght nnz


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
