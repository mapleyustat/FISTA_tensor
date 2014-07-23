function [opt, objV] = proxF_nuc(x,tau,msize) 

% reshape x into matrix form first.
% At last return a vector

x        =  reshape(x,msize);
[u s v]  =  svd(x,'econ')   ;
[s,objV] =  proxF_l1(s,tau) ;
opt      =  u*s*v'          ;
opt      =  opt(:)          ;