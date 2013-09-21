function [data] = lasso()
    randn('seed',0); rand('seed',0);
    m = 20;%2500;
    n = 100;%50000;
    A = randn(m,n); b = randn(m,1);
    gamma = 1;
%{
    tic
    cvx_solver sedumi
    cvx_begin
        variable x(n)
        minimize (sum_square(A*x - b) + gamma * norm(x,1))
    cvx_end
    sedumi_time = toc

    save sedumi_time sedumi_time m n
%}

    dims = struct('m',m,'n',n);
    params = struct('A',A,'b',b,'gamma',gamma);
    
    socp_data = prob_to_socp(params, dims);
    
%    ecos(data.c, data.G, data.h, data.dims);
    
    mb = nnz(socp_data.G)*8/(1024*1024);
    fprintf('total data is %f mb\n',mb);

    data.A = socp_data.G;
    data.b = socp_data.h;
    data.c = socp_data.c;
    cones = socp_data.dims;
    cones.f = 0;

    socp_data

    params.VERBOSE = 1;
    write_pdos_data(data, cones, params, 'lasso_data2');

end

function [data] = prob_to_socp(params, dims)
    % PROB2SOCP: maps PARAMS into a struct of SOCP matrices
    % Where input struct PARAMS has the following fields:
    %   'A' has shape Matrix(m,n)
    %   'b' has shape Vector(m)
    %   'gamma' has shape Scalar()

    p = 0; m = (4 + dims.m + 2*dims.n); n = (2 + 2*dims.n);
    c = zeros(n,1);
    h = zeros(m,1);
    b = zeros(p,1);
    Gi = []; Gj = []; Gv = [];
    Ai = []; Aj = []; Av = [];

    cones.l = 2*dims.n;
    cones.q = [(1 + dims.m); 3];
    cones.s = [];

    % stuffing the objective vector
    c((2 + 2*dims.n):(2 + 2*dims.n)) = 1;
    c((2 + dims.n):(1 + 2*dims.n)) = params.gamma * ones(dims.n,1);

    % for the constraint -1*_t2 + -1*x <= 0
    Gi = [Gi; (0:(-1 + dims.n))'];
    Gj = [Gj; (1:dims.n)'];
    Gv = [Gv; -1*ones(dims.n,1)];
    Gi = [Gi; (0:(-1 + dims.n))'];
    Gj = [Gj; ((1 + dims.n):2*dims.n)'];
    Gv = [Gv; -1*ones(dims.n,1)];

    % for the constraint x + -1*_t2 <= 0
    Gi = [Gi; (dims.n:(-1 + 2*dims.n))'];
    Gj = [Gj; (1:dims.n)'];
    Gv = [Gv; 1*ones(dims.n,1)];
    Gi = [Gi; (dims.n:(-1 + 2*dims.n))'];
    Gj = [Gj; ((1 + dims.n):2*dims.n)'];
    Gv = [Gv; -1*ones(dims.n,1)];

    % for the SOC constraint norm([A*x + -1*b]) <= _t0
    h((2 + 2*dims.n):(1 + dims.m + 2*dims.n)) = -(params.b);
    Gi = [Gi; (1 + 2*dims.n) + mod(find(params.A)-1,size(params.A,1))];
    Gj = [Gj; 1 + floor((find(params.A)-1)/size(params.A,1))];
    Gv = [Gv; nonzeros(-params.A)];
    Gi = [Gi; 2*dims.n];
    Gj = [Gj; 0];
    Gv = [Gv; -1];

    % for the SOC product constraint norm(1 + -1*_t1, 2.0*_t0) <= 1 + _t1
    Gi = [Gi; (3 + dims.m + 2*dims.n)];
    Gj = [Gj; 0];
    Gv = [Gv; -2.0];
    h((3 + dims.m + 2*dims.n):3:(4 + dims.m + 2*dims.n)) = 1;
    Gi = [Gi; (2 + dims.m + 2*dims.n)];
    Gj = [Gj; (1 + 2*dims.n)];
    Gv = [Gv; 1];
    h((2 + dims.m + 2*dims.n):3:(4 + dims.m + 2*dims.n)) = 1;
    Gi = [Gi; (1 + dims.m + 2*dims.n)];
    Gj = [Gj; (1 + 2*dims.n)];
    Gv = [Gv; -1];

    % Convert from sparse triplet to column compressed format.
    % Also convert from 0 indexed to 1 indexed.
    A = sparse(Ai+1, Aj+1, Av, p, n);
    G = sparse(Gi+1, Gj+1, Gv, m, n);

    % Build output
    data = struct('c', c, 'b', b, 'h', h, 'G', G, 'A', A, 'dims', cones);
end
