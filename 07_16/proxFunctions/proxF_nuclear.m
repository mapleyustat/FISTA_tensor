function [x, objV] = proxF_nuclear(x,t) 

[u s v] = svd(x,'econ');

[s,objV] = proxF_l1(s,t);

x = u*s*v';
