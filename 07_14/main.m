clear all
close all
clc

n1 = 5;
n2 = 5;
n3 = 3;
p  = 0.5;
lambda = 1e-0;

m = round(p*n1*n2*n3) ;
x = randn(n1,n2,n3)   ;
[U S V] = ntsvd(x)    ;
thr     = S(3,3,1)    ;
S       = S.*(S>=thr) ;
x  = tprod(tprod(U,S),tran(V));
A = randn(m,n1*n2*n3) ;
b = A*x(:)            ;
QUIET = 0;

tensor_size = [n1 n2 n3] ;
% x_opt = FISTA_const( A , b , [n1 n2 n3] , lambda , QUIET);
% overall parameters
tol    =  1e-1                          ;

[~,n]  =  size(A)                       ;
AA     =  A'*A                          ;
L      =  2*max(eig(AA))                ;
B      =  ones(size(AA)) - 2*AA/L       ;
C      =  2/L*A'*b                      ;

% step 0
x      =  5*randn(n,1)                    ; 
y      =  x                             ;
t      =  1                             ;
iter   =  0                             ;
max_ite=  10                             ; 

% step 1
if ~QUIET
        fprintf('%3s\t%7s\t%16s\t%7s\n',...
            'iter','t','TNN','gap')         ;
end
    
while norms(A*x-b) > tol && iter < max_ite
    iter   =  iter+1                    ;
    x_old  =  x                         ;
    t_old  =  t                         ;
    
    Tm     =  B*y+C                     ;
    T      =  reshape(Tm,tensor_size)   ;
    Tf     =  fft(T,[],3)               ;
    Tfbm   =  t_blkdiag(Tf)             ;
    
    [x TNN]=  pL( Tfbm , lambda/L)      ;
    xf     =  t_iblkdiag(x,tensor_size)    ;
    xt     =  ifft(xf,[],3)              ;
    x      =  real(xt(:))                   ;
    
    t      =  (1+sqrt(1+4*t^2))/2       ;
    y      =  x + (t_old-1)/t*(x-x_old) ;
    gap    =  norms(A*x-b)              ;
    
    if ~QUIET
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\n',...
            iter,t,TNN,gap)               ;
    end

end

x_opt  = reshape(x,tensor_size)          ;