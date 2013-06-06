function [x,s,y,status, result] = pdos(data,cone,params)
% cone solver, solves:
%
% min. c'x
% subject to Ax + s = b
% s \in K
%
% where x \in R^n, s \in R^m
%
% K is product of cones in this particular order:
% zero cone, lp cone, second order cone(s)
%
% data must consist of data.A, data.b, data.c, where A,b,c used as above
%
% cone struct must consist of:
% cone.f, length of zero cone (for equality constraints)
% cone.lp, length of lp cone
% cone.k_soc, number of second-order cones (SOC)
% cone.q, array of SOC lengths
% (thus m = cones.zero+cone.lp+sum(cone.q) = size(A,1))
%
% cone struct is only used in proj_cone, to add other cones
% simply add the relevant size data to the cone struct and edit the
% proj_cone method to include projection onto the new cone AND
% the dual cone (current implementation is second-order cones only)

% params struct consists of the following fields
% (here set to default settings):
MAX_ITERS = 2000; % maximum num iterations for admm
EPS_ABS   = 1e-3; % quitting tolerances
VERBOSE   = true; % verbosity level

% conjugate gradient (CG) settings:
USE_CG       = false;  % whether to use CG
CG_MAX_ITERS = 20; % max iterations for CG
CG_TOL       = 1e-4; % max CG quitting tolerance
CG_VERBOSE   = false; % CG prints summary

NORMALIZE = 0;
ALPHA = 1.0;
WARMSTART = false;

%%
if nargin==3
    if isfield(params,'MAX_ITERS');MAX_ITERS = params.MAX_ITERS;end
    if isfield(params,'EPS_ABS');EPS_ABS = params.EPS_ABS;end    
    if isfield(params,'WARMSTART'); WARMSTART = params.WARMSTART; end
    
    % CG:
    if isfield(params,'USE_CG'); USE_CG = params.USE_CG; end
    if isfield(params,'CG_MAX_ITS');CG_MAX_ITERS = params.CG_MAX_ITS;end
    if isfield(params,'CG_TOL');CG_TOL = params.CG_TOL;end
    if isfield(params,'CG_VERBOSE');CG_VERBOSE = params.CG_VERBOSE;end
    
    % VERBOSITY
    if isfield(params,'VERBOSE'); VERBOSE = params.VERBOSE; end
    % DLIZE
    if isfield(params,'DLIZE'); NORMALIZE = params.NORMALIZE; end
    
    if isfield(params,'ALPHA'); ALPHA = params.ALPHA; end
end


%%
n = length(data.c);
m = length(data.b);

if (CG_MAX_ITERS > m)
    CG_MAX_ITERS = m;
end

%%
D = ones(m,1);  % D
E = ones(n,1);  % E
ratio = 1.0;    % 1 / sqrt( mu )

% store these before Dlizing
bnorm = norm(data.b,'inf');
cnorm = norm(data.c,'inf');


if NORMALIZE
    ratio = 1000;
    N = 3;

    p = cone.f + cone.l + length(cone.q);

    d1 = ones(p+1,1);
    e1 = ones(n+1,1);
    e = e1;

% we implement the following lines of code using sparse matrix manipulation
% in what follows. when finished, As = Atilde
%{
    Atilde = sparse(p+1,n+1);
    Atilde(1:cone.f+cone.l,end) = sparse(data.b(1:cone.f+cone.l));
    Atilde(end,1:end-1) = sparse(data.c');
    Atilde(1:cone.f+cone.l,1:end-1) = data.A(1:cone.f+cone.l,:);
    
    ind = cone.f + cone.l;
    for j = 1:cone.k_soc,
        soc_sz = cone.q(j);
        Atilde(cone.f+cone.l+j,:) = [norms(data.A(ind+1:ind+soc_sz,:)) norm(data.b(ind+1:ind+soc_sz))] ;
        ind = ind + soc_sz;
    end
%}
    
    [AI,AJ,AV] = find(data.A);
    ind = AI <= cone.f+cone.l;
    J = AJ(ind);
    V = AV(ind);
    I = AI(ind);
    
    [bI,bJ,bV] = find(data.b);
    ind = bI <= cone.f+cone.l;

    I = [I; bI(ind)];
    J = [J; bJ(ind)+n];
    V = [V; bV(ind)];
    
    ind = cone.f + cone.l;
    for j = 1:cone.k_soc,
        soc_sz = cone.q(j);
        tmp = sparse( [norms(data.A(ind+1:ind+soc_sz,:), inf) norm(data.b(ind+1:ind+soc_sz), inf)] );%/sqrt(soc_sz);
        [tI,tJ,tV] = find(tmp);
        I = [I; (cone.f+cone.l+j)*ones(length(tI),1)];
        J = [ J; tJ' ];
        V = [ V; tV' ];
        ind = ind + soc_sz;
    end
    
    [cI,cJ,cV] = find(data.c);
    I = [I; cJ+p];
    J = [J; cI];
    V = [V; cV];
    
    H = sparse(I,J,V, p+1,n+1);

    for i = 1:N,
        H = H*spdiags(e,0,n+1,n+1);        
        d = sqrt( 1./ norms(As',inf)' );
        
        H = spdiags(d,0,p+1,p+1)*H;
        e = sqrt( 1./norms(As,inf)' ) ;

        d1 = d1.*d;
        e1 = e1.*e;
    end
    
    % extend the blocks
    dd = zeros(n,1);
    dd(1:cone.f+cone.l) = d1(1:cone.f+cone.l);
    ind = cone.f + cone.l;
    for j = 1:length(cone.q),
        soc_sz = cone.q(j);
        dd(ind+1:ind+soc_sz,1) = d1(cone.f+cone.l+j);
        ind = ind + soc_sz;
    end
    
    D = dd;
    E = e1(1:end-1);
    
    
    data.b = D.*data.b;
    data.c = E.*data.c;
    
    lambda = norm(data.b)/norm(data.c)   
    
    E = E * ratio;    % assumes mu chosen so mu/lambda = 1e-6 (so that it scales with mu)
    data.c = data.c * ratio;
    data.A = spdiags(D,0,m,m)*data.A*spdiags(E,0,n,n);

else

    lambda = norm(data.b)/norm(data.c)     

end

W=sparse([speye(m,m) data.A;
 data.A' -speye(n,n)]);
q = amd(W);
[L,ldl_diag] = ldlsparse(W,q);
LD = L + ldl_diag;

if VERBOSE
    disp('Solving cone program')
    tic
end

headings = {('Iter'), ('ni(Ax + s - b)'), ('ni(A''y + c)'), ...
    ('   c''x   '), ('  -b''y   '), ('   eta   ')};
if VERBOSE
    fprintf('\n');
    line_len = 0;
    for i = 1:length(headings)-1
        fprintf('%*s | ', length(headings{i}), headings{i});
        line_len = line_len + length(headings{i}) +3;
    end
    fprintf('%*s\n', length(headings{end}), headings{end});
    line_len = line_len + length(headings{end});

    line = repmat('-',1,line_len);
    fprintf('%s\n',line);
end

status = 'Max Iters';


x = zeros(n,1); s = zeros(m,1); y = zeros(m,1);

if WARMSTART
    if isfield(params,'X_INIT'); x = params.X_INIT ./ E; end
    if isfield(params,'Y_INIT'); y = params.Y_INIT ./ D; end
    if isfield(params,'S_INIT'); s = params.S_INIT .* D; end
end

result.primal_residual = zeros(MAX_ITERS,1);
result.dual_residual = zeros(MAX_ITERS,1);
result.gap = zeros(MAX_ITERS,1);
result.eps_pri = zeros(MAX_ITERS,1);
result.eps_dual = zeros(MAX_ITERS,1);
result.eps_gap = zeros(MAX_ITERS,1);
result.stop = 0;

result.total_cg = 0;
for i=1:MAX_ITERS
    
    %% ADMM iterates
    x_old = x;
    
    % project onto linear equality constraints
    if USE_CG
        % CG_TOL is multiplied by ratio, since error in "x" space is
        % amplified by "ratio"
        [x,iter] = pcg_custom(data.A, data.b, x_old, x_old-lambda*data.c, s + lambda*y, CG_MAX_ITERS, CG_TOL*ratio, CG_VERBOSE);
        result.total_cg = result.total_cg + iter;
    else
        rhs = [data.b-lambda*y-s; (lambda*data.c - x_old)];
        x = ldlsolve(LD,rhs(q));x(q)=x;
        % tmp = x(1:m);
        x = x(m+1:end);
    end

    % overreleax
    s_over = ALPHA*(data.b - data.A*x) + (1-ALPHA)*s;
    
    % primal cone projection
    s = proj_cone(s_over - lambda*y,cone);
    
    % dual variable update
    y = y + (1.0/lambda)*(s - s_over);
    
    %% stopping conditions:
    p_res = norm( (data.A*x + s - data.b)./D );
    d_res = norm( (data.A'*y + data.c)./E);
    eta = norm(data.c'*x + data.b'*y);
    
    p_obj = data.c'*x;
    d_obj = -data.b'*y;
    
    result.primal_residual(i) = p_res;
    result.dual_residual(i) = d_res;
    result.gap(i) = eta;
    result.eps_pri(i) = EPS_ABS*(1 + bnorm);
    result.eps_dual(i) = EPS_ABS*(1 + cnorm);
    result.eps_gap(i) = EPS_ABS*(1 + abs(p_obj) + abs(d_obj));
    
    if (p_res < result.eps_pri(i)) && (d_res < result.eps_dual(i)) && (eta < result.eps_gap(i))
        status = 'Solved';
        if result.stop == 0
            result.stop = i;
        end
        break
    end
    
    
    if VERBOSE && mod(i-1,10)==0
        print_summary(headings,i,p_res,d_res,p_obj,d_obj,eta);
    end
end

if VERBOSE
    print_summary(headings,i,p_res,d_res,p_obj,d_obj,eta);
end


if VERBOSE
    fprintf('Took %i iterations\n',i)
    toc
end

x = E.*x;
s = s ./ D;
y = D .* y;

end

function print_summary(headings,iter,p_res,d_res,p_obj,d_obj,eta)
    fprintf('%*i | ', length(headings{1}), iter-1);
    fprintf('%*.2e   ', length(headings{2}), full(p_res));
    fprintf('%*.2e   ', length(headings{3}), full(d_res));
    fprintf('%*.2e   ', length(headings{4}), full(p_obj));
    fprintf('%*.2e   ', length(headings{5}), full(d_obj));
    fprintf('%*.2e\n', length(headings{6}), full(eta));
end

function z = proj_cone(z,c)
    free_len = c.f;
    lp_len = c.l;
    k_soc = length(c.q);
    ns_soc = c.q;
    % zero cone, (dual = free cone)
    z(1:free_len) = zeros(free_len,1);  
    
    % lp cone
    z(free_len+1:lp_len+free_len) = pos(z(free_len+1:lp_len+free_len));
    
    % SOCs
    idx=lp_len+free_len;
    idxs = [idx; idx + cumsum(ns_soc)];
    
    for i=1:k_soc
        z(idxs(i)+1:idxs(i)+ns_soc(i)) = proj_soc(z(idxs(i)+1:idxs(i)+ns_soc(i)));
    end
end

function z = proj_soc(tt)
    v1=tt(1);v2=tt(2:end);
    if norm(v2)<=-v1
        v2=zeros(length(v2),1);v1=0;
    elseif norm(v2)> abs(v1)
        v2=0.5*(1+v1/norm(v2))*v2;
        v1=norm(v2);
    end
    z=[v1;v2];
end

function [x,i] = pcg_custom(A,b,x,u,v,MAX_ITS,TOL,VERBOSE)
% diag preconditioner
M = 1;

NEWTOL = TOL*TOL;
At = A';
r = u - x - At*(A*x - b + v);
z = M\r;
p = z;
rsold=z'*r;
if (r'*r) >= NEWTOL

    for i=1:MAX_ITS
        iAp=p + At*A*p;
        alpha=rsold/(p'*iAp);
        x=x+alpha*p;
        r=r-alpha*iAp;
        if (r'*r) < NEWTOL*1e7
            if VERBOSE
                fprintf('CG took %i iterations to reach tolerance %4f\n',i,TOL)
            end
            return;
        end
        z = M\r;
        rsnew=z'*r;

        p=z+(rsnew/rsold)*p;
        rsold=rsnew;
    end
    
    if VERBOSE
        fprintf('CG did not converge within %i iterations\n',MAX_ITS)
    end
end
end
