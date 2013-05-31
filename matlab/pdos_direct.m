function [x, s, y, status] = pdos_direct(data, cone, params)			    
%PDOS_DIRECT Operator-splitting method for solving cone problems (direct)
%
% This implements a cone solver. It solves:
%
% min. c'x
% subject to Ax + s = b
% s \in K
%
% where x \in R^n, s \in R^m
%
% K is product of cones in this particular order:
% free cone, lp cone, second order cone(s)
%
% data must consist of data.A, data.b, data.c, where A,b,c used as above.
% You must take care to ensure that data.A is sparse while data.b and data.c
% are dense. 
%  
% cone struct must consist of:
% cone.f, length of free cone (for equality constraints)
% cone.l, length of lp cone
% cone.q, array of SOC lengths
%
% cone struct is only used in proj_cone, to add other cones
% simply add the relevant size data to the cone struct and edit the
% proj_cone method to include projection onto the new cone AND
% the dual cone (current implementation is second-order cones only)
%
% Necessary fields in the params struct are:
%   ALPHA       : over-relaxation parameter, between (0,2).
%   MAX_ITERS   : maximum number of ADMM iterations.
%   EPS_ABS     : accuracy of solution
%   CG_MAX_ITS  : maximum number of CG iterations
%   CG_TOL      : tolerance of CG
%   VERBOSE     : verbosity level (0 or 1)
error ('pdos_direct mexFunction not found') ;
