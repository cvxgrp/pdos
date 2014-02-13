function shim = cvx_pdos( shim )

global cvx___ %#ok
if ~isempty( shim.solve ),
    return
end
if isempty( shim.name ),
    [ fs, ps, int_path, mext ] = cvx_version; %#ok
    % DIRECT version
    fname = [ 'pdos_direct.', mext ];
    % INDIRECT version (uncomment to use)
    % fname = [ 'pdos_indirect.', mext ];
    shim.name = 'pdos';
    shim.dualize = true;
    flen = length(fname);
    fpaths = which( fname, '-all' );
    if ~iscell(fpaths),
      fpaths = { fpaths };
    end
    old_dir = pwd;
    oshim = shim;
    shim = [];
    for k = 1 : length(fpaths),
        fpath = fpaths{k};
        if ~exist( fpath, 'file' ) || any( strcmp( fpath, fpaths(1:k-1) ) ),
            continue
        end
        new_dir = fpath(1:end-flen-1);
        cd( new_dir );
        tshim = oshim;
        tshim.fullpath = fpath;
        tshim.version = 'unknown';
        tshim.location = new_dir;
        if isempty( tshim.error ),
            tshim.check = @check;
            tshim.solve = @solve;
            tshim.eargs = {};
            if k ~= 1,
                tshim.path = [ new_dir, ps ];
            end
        end
        shim = [ shim, tshim ]; %#ok
    end
    cd( old_dir );
    if isempty( shim ),
        shim = oshim;
        shim.error = 'Could not find a PDOS installation.';
    end
else
    shim.check = @check;
    shim.solve = @solve;
end

function found_bad = check( nonls, sdp, mfunc ) %#ok
found_bad = false;
if ~sdp,
    for k = 1 : length( nonls ),
        if any( strcmp( nonls(k).type, { 'semidefinite', 'hermitian-semidefinite' } ) ) && size(nonls(k).indices,1) > 4
            warning( 'CVX:PDOS:Semidefinite', 'This nonlinearity requires use of semidefinite cones which PDOS does not support.\n%s', ...
                'You will need to use a different solver for this model.' );
            found_bad = true;
        end
    end
end

function [ x, status, tol, iters, y, z ] = solve( At, b, c, nonls, quiet, prec, settings )

n = length( c );
m = length( b );
K = struct( 'f', 0, 'l', 0, 'q', [], 'r', [] );
reord = struct( 'n', 0, 'r', [], 'c', [], 'v', [] );
reord = struct( 'f', reord, 'l', reord, 'a', reord, 'q', reord, 'r', reord );
reord.f.n = n;
zinv = [];
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
            zinv = [ zinv ; temp(:) ]; %#ok
        else
            temp = temp( [ end, 1 : end - 1 ], : );
            reord.q.r = [ reord.q.r ; temp(:) ];
            reord.q.c = [ reord.q.c ; reord.q.n + ( 1 : nnv )' ];
            reord.q.v = [ reord.q.v ; ones(nnv,1) ];
            reord.q.n = reord.q.n + nnv;
            K.q = [ K.q, nn * ones( 1, nv ) ];
        end
    elseif isequal( tt, 'semidefinite' ),
        if nn == 3,
            temp = temp( [1,1,3,3,2], : );
            tempv = [1;1;1;-1;1] * ones(1,nv);
            tempc = bsxfun(@plus,[1;2;1;2;3],3*(0:nv-1));
            reord.q.r = [ reord.r.r ; temp(:) ];
            reord.q.c = [ reord.r.c ; reord.r.n + tempc(:) ];
            reord.q.v = [ reord.r.v ; tempv(:) ];
            reord.q.n = reord.r.n + nnv;
            K.q = [ K.q, 3 * ones( 1, nv ) ];
            temp = temp([1,3],:);
            zinv = [ zinv ; temp(:) ]; %#ok
        else
            error( 'CVX:SolverIncompatible', 'PDOS does not support semidefinite variables larger than 2x2.' );
        end
    elseif isequal( tt, 'hermitian-semidefinite' ),
        if nn == 4,
            temp = temp( [1,1,4,4,2,3], : );
            tempv = [1;1;1;-1;1;1] * ones(1,nv);
            tempc = bsxfun(@plus,[1;2;1;2;3;4],4*(0:nv-1));
            reord.q.r = [ reord.r.r ; temp(:) ];
            reord.q.c = [ reord.r.c ; reord.r.n + tempc(:) ];
            reord.q.v = [ reord.r.v ; tempv(:) ];
            reord.q.n = reord.r.n + nnv;
            K.q = [ K.r, q * ones( 1, nv ) ];
            temp = temp([1,3],:);
            zinv = [ zinv ; temp(:) ]; %#ok
        else
            error( 'CVX:SolverIncompatible', 'PDOS does not support semidefinite variables larger than 2x2.' );
        end
    else
        error( 'Unsupported nonlinearity: %s', tt );
    end
end
if reord.f.n > 0,
    reord.f.r = ( 1 : n )';
    reord.f.r( [ reord.l.r ; reord.a.r ; reord.q.r ] ) = [];
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
reord = sparse( ...
    [ reord.f.r ; reord.l.r ; reord.a.r ; reord.q.r ], ...
    [ reord.f.c ; reord.l.c ; reord.a.c ; reord.q.c ], ...
    [ reord.f.v ; reord.l.v ; reord.a.v ; reord.q.v ], ...
    n, n_out );

At = reord' * At;
c  = reord' * c;

pars.ALPHA = 1.8;
pars.MAX_ITERS = 10000;
pars.EPS_ABS = 1e-3;
if quiet
    pars.VERBOSE = 0;
else
    pars.VERBOSE = 1;
end
pars.NORMALIZE = 1;

pars.CG_MAX_ITS = 20;
pars.CG_TOL = 1e-4;

% add a row if the problem is unconstrained
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

if ~isempty(zinv),
    z(zinv) = z(zinv) * 0.5;
end
