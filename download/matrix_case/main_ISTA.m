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
n1     = 50 ;
rank   = 5 ;
p      = 0.7;

X      =  randn(n1,rank)*randn(rank,n1) ;
msize  =  size(X)                        ;
A      =  randn(round(p*n1^2),n1^2)     ;
b      =  A*X(:)                        ;

% overall parameters
tol    =  1e-1                          ;
lambda =  10                          ;
QUIET  =  0                             ;

L      =  0.00001*rand                          ;
eta    =  1.5                       ;

[~,n]  =  size(A)                       ;
AA     =  A'*A                          ;
L      =  2*max(eig(AA))                ;
B      =  ones(size(AA)) - 2*AA/L       ;
C      =  2*A'*b/L                      ;


% step 0
x      =  5*randn(n,1)                  ;
xm     = reshape(x,msize) ;
iter   =  0                             ;
max_ite=  100                            ; 
err    =  10                            ;

imax   = 1000;

f = @(x) norm(A*x(:)-b)^2 ;
g = @(x) lambda*sum(svd(x)) ;

% step 1
if ~QUIET
        fprintf('%3s\t%16s\t%7s\n',...
            'iter','TNN','err')         ;
end
    
while err > tol && iter < max_ite
    
    for i = 1:imax
        disp(i)
        L = L*eta^i;
        
        % compute F(pL(x_{k-1}))
        vtmp = B*x+C;
        Mtmp = reshape(vtmp,msize);
        pLx  =  pL(Mtmp,2*lambda/L);
        Left = f(pLx) + g(pLx)     ;
        
        
        % compute Q_L
        Right = f(xm)+(pLx(:)'-x')*(2*AA*x-2*A'*b)+L/2*norm(x-pLx(:))^2 + g(pLx);
        
        if Left <= Right
            break;
        end
    end
    
    
    
    iter   =  iter+1                    ;
    
    ve      =  B*x+C                    ;
    Ma      =  reshape(ve,msize)        ;
    s       =  svd(Ma)                  ;
    [xm TNN]=  pL( Ma , 2*lambda/L)     ;
    err     =  norms(A*xm(:)-b)         ;

%     if ~QUIET
%         fprintf('%3d\t%10.4f\t%10.4f\n',...
%             iter,TNN,err)               ;
%     end

end

x_opt  = xm ; 