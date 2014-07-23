% ****Solve the problem using FISTA with constant L****
%
%                     min ||X||_*
%             subject to AX(:) = b
%  ==> 
%                 min ||AX(:)-b||+lambda||X||_*
%
% *****************************************************
% by Jamie Zemin Zhang
% 07/15/2014
%

clear all
close all
clc


% problem set up
n1     = 20 ;
rank   = 4 ;
p      = 0.4;

X      =  randn(n1,rank)*randn(rank,n1) ;
msize  = size(X)                        ;
A      =  randn(round(p*n1^2),n1^2)     ;
b      =  A*X(:)                        ;

% overall parameters
tol    =  1e-1                          ;
lambda =  1e-2                          ;
QUIET  =  0                             ;


[~,n]  =  size(A)                       ;
AA     =  A'*A                          ;
L      =  2*max(eig(AA))                ;
B      =  ones(size(AA)) - 2*AA/L       ;
C      =  2*A'*b/L                      ;


% step 0
x      =  5*randn(n,1)                  ; 
y      =  x                             ;
t      =  1                             ;
iter   =  0                             ;
max_ite=  10                            ; 
err    =  10                            ;

% step 1
if ~QUIET
        fprintf('%3s\t%7s\t%16s\t%7s\n',...
            'iter','t','TNN','err')         ;
end
    
while err > tol && iter < max_ite
    iter   =  iter+1                    ;
    x_old  =  x                         ;
    t_old  =  t                         ;
    
    ve      =  B*y+C                     ;
    Ma      =  reshape(ve,msize)   ;
    
    [xm TNN]  =  pL( Ma , 2*lambda/L)      ;
    
    x       =  xm(:)                     ;
    t       =  (1+sqrt(1+4*t^2))/2       ;
    y       =  x + (t_old-1)/t*(x-x_old) ;
    err     =  norms(A*x-b)              ;
    
    if ~QUIET
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\n',...
            iter,t,TNN,err)               ;
    end

end

x_opt  = xm ; 