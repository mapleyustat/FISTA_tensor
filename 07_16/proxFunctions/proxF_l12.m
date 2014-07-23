function [x, objV] = proxF_l12(x,t) 
%% Usage proxF_l1(x,t)
% ProxF_l12(x,t) applies L1 shrinkage to 2-dimensional matrix x for
% a shrinkage factor of t. Shrinkage is applied along columns of x.
%Returns the shrunk object x.
%
%
% INPUTS:
%
% x: 2-dimensional matrix
% t: shrinkage threshold
%
% OUTPUTS:
%
% x: n-dimensional matrix after shrinkage.
%
v = sqrt( sum(x.^2,2) );

s = 1 - 1 ./ max( v ./ ( t ), 1 );
m = length(s);
x = spdiags(s,0,m,m)*x;

objV = sum(sqrt( sum(x.^2,2) ));