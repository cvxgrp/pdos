function shim = cvx_pdos( shim )

% CVX_SOLVER_SHIM	pdos interface for CVX.
%   This procedure returns a 'shim': a structure containing the necessary
%   information CVX needs to use this solver in its modeling framework.

if ~isempty( shim.solve ),
    return
end
if isempty( shim.name ),
    [ fs, ps ] = cvx_version; %#ok
    temp = strfind( shim.spath, fs );
    shim.name    = 'pdos';
    shim.dualize = true;   % don't allow dualization, yet
    shim.path    = [ shim.spath(1:temp(end-1)), 'pdos', ps ];
end
if isempty( shim.error ),
    shim.check = @check;
    shim.solve = @solve;
else
    shim.check = [];
    shim.solve = [];
end

function found_bad = check( nonls ) %#ok
found_bad = false;

function [ x, status, tol, iters, y, z ] = solve( At, b, c, nonls, quiet, prec, settings )

n = length( c );
m = length( b );
K = struct( 'f', 0, 'l', 0, 'q', [], 's', [], 'scomplex', [], 'ycomplex', [] );
reord = struct( 'n', 0, 'r', [], 'c', [], 'v', [] );
reord = struct( 'f', reord, 'l', reord, 'a', reord, 'q', reord, 's', reord, 'h', reord );
reord.f.n = n;
for k = 1 : length( nonls ),
    temp = nonls( k ).indices;
    nn = size( temp, 1 );
    nv = size( temp, 2 );
    nnv = nn * nv;
    tt = nonls( k ).type;
    reord.f.n = reord.f.n - nnv;
    if strncmp( tt, 'i_', 2 ),
        error( 'PDOS does not support integer variables.' );
    elseif nn == 1 || isequal( tt, 'nonnegative' ),
        reord.l.r = [ reord.l.r ; temp(:) ];
        reord.l.c = [ reord.l.c ; reord.l.n + ( 1 : nnv )' ];
        reord.l.v = [ reord.l.v ; ones( nnv, 1 ) ];
        reord.l.n = reord.l.n + nnv;
    elseif isequal( tt, 'lorentz' ),
        if nn == 2,
            rr = [ temp ; temp ];
            cc = reshape( floor( 1 : 0.5 : 2 * nv + 0.5 ), 4, nv );
            vv = [1;1;-1;1]; vv = vv(:,ones(1,nv));
            reord.a.r = [ reord.a.r ; rr(:) ];
            reord.a.c = [ reord.a.c ; cc(:) + reord.a.n ];
            reord.a.v = [ reord.a.v ; vv(:) ];
            reord.a.n = reord.a.n + nnv;
        else
            temp = temp( [ end, 1 : end - 1 ], : );
            reord.q.r = [ reord.q.r ; temp(:) ];
            reord.q.c = [ reord.q.c ; reord.q.n + ( 1 : nnv )' ];
            reord.q.v = [ reord.q.v ; ones(nnv,1) ];
            reord.q.n = reord.q.n + nnv;
            K.q = [ K.q, nn * ones( 1, nv ) ];
        end
    elseif isequal( tt, 'semidefinite' ),
        % need to convert 2x2 into SOCs
        error( 'PDOS does not support semidefinite programming.' );
    elseif isequal( tt, 'hermitian-semidefinite' ),
        error( 'PDOS does not support semidefinite programming.' );
    else
        error( 'Unsupported nonlinearity: %s', tt );
    end
end
if reord.f.n > 0,
    reord.f.r = ( 1 : n )';
    reord.f.r( [ reord.l.r ; reord.a.r ; reord.q.r ; reord.s.r ] ) = [];
    reord.f.c = ( 1 : reord.f.n )';
    reord.f.v = ones(reord.f.n,1);
end
n_d = max( m - n - reord.f.n + 1, isempty( At ) );
if n_d,
    reord.l.n = reord.l.n + n_d;
end
K.f = reord.f.n;
K.l = reord.l.n + reord.a.n;
n_out = reord.f.n;
reord.l.c = reord.l.c + n_out; n_out = n_out + reord.l.n;
reord.a.c = reord.a.c + n_out; n_out = n_out + reord.a.n;
reord.q.c = reord.q.c + n_out; n_out = n_out + reord.q.n;
reord.s.c = reord.s.c + n_out; n_out = n_out + reord.s.n;
reord = sparse( ...
    [ reord.f.r ; reord.l.r ; reord.a.r ; reord.q.r ; reord.s.r ], ...
    [ reord.f.c ; reord.l.c ; reord.a.c ; reord.q.c ; reord.s.c ], ...
    [ reord.f.v ; reord.l.v ; reord.a.v ; reord.q.v ; reord.s.v ], ...
    n, n_out );

At = reord' * At;
c  = reord' * c;
%pars.free = K.f > 1 && nnz( K.q );
%pars.eps     = prec(1);
%pars.bigeps  = prec(3);
%if quiet,
%    pars.fid = 0;
%end

pars.ALPHA = 1.8;
pars.MAX_ITERS = 10000;
pars.EPS_ABS = 5e-3;
if quiet
    pars.VERBOSE = 0;
else
    pars.VERBOSE = 1;
end
pars.NORMALIZE = 1;

pars.CG_MAX_ITS = 20;
pars.CG_TOL = 1e-4;

add_row = isempty( At );
if add_row,
    K.f = K.f + 1;
    At = sparse( 1, 1, 1, n_out + 1, 1 );
    b = 1;
    c = [ 0 ; c ];
end

data = struct('A', At, 'b', full(c), 'c', -full(b));

% transpose the cone if needed, since pdos_direct assumes cone is a column
% vector argument
[m1,n1] = size(K.q);
if m1 == 1,
    K.q = K.q';
end

% coneOS solves the dual formulation of sedumi canonical form
% IMPORTANT!!!: At appears to be A', sedumi seems not to care...
% USE DIRECT METHOD:
[ yy, ss, xx, status ] = cvx_run_solver( @pdos_direct, data,K, pars, 'xx', 'ss', 'yy', 'info', settings, 4 );

% UNCOMMENT TO USE CG:
%[ yy, xx, status ] = cvx_run_solver( @pdos_indirect, At, full(c),-full(b),K, pars, 'xx', 'yy', 'info', settings, 4 );


if add_row,
    xx = xx(2:end);
    yy = zeros(0,1);
    At = zeros(n_out,0);
    % b  = zeros(0,1);
    c  = c(2:end);
end
tol = c'*xx + b'*yy;
iters = 0;
xx = full( xx );
yy = full( yy );

if strcmp(status,'Max Iters');
    status = 'Indeterminate';
    x = NaN * ones( n, 1 );
    y = NaN * ones( m, 1 );
    z = NaN * ones( n, 1 );
    if add_row, y = zeros( 0, 1 ); end
else
    x = real( reord * xx );
    y = yy;
    z = real( reord * ( c - At * yy ) );
    if add_row, y = zeros( 0, 1 ); end
end
