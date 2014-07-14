function x_opt = FISTA_const( A , b , size , QUIET)
%
% ****Solve the problem using FISTA with constant L****
%
%                     min ||X||_{TNN}
%             subject to AX(:) = b
%  ==> 
%                 min ||AX(:)-b||+lambda||X||_{TNN}
%
% *****************************************************
% by Jamie Zemin Zhang
% 07/14/2014
%

% overall parameters
lambda =  1e-3                          ;
tol    =  1e-2                          ;

[~,n]  =  size(A)                       ;
AA     =  A'*A                          ;
L      =  2*max(eig(AA))                ;
B      =  ones(size(AA)) - 2*AA/L       ;
C      =  2/L*A'*b                      ;

% step 0
x      =  randn(size(A,2),1)                    ; 
y      =  x                             ;
t      =  1                             ;
iter   =  0                             ;

% step 1
if ~QUIET
        fprintf('%3s\t%10s\t%10s\n',...
            'iter','TNN','gap')         ;
end
    
while norms(A*x-b) > tol
    iter   =  iter+1                    ;
    x_old  =  x                         ;
    t_old  =  t                         ;
    T      =  B*y+C                     ;
    Tf     =  fft(T,[],3)               ;
    Tfbm   =  t_blkdiag(Tf)             ;
    
    [x TNN]=  pL( Tfbm , lambda/L)      ;
    t      =  (1+sqrt(1+4*t^2))/2       ;
    y      =  x + (t_old-1)/t*(x-x_old) ;
    gap    =  norms(A*x-b)              ;
    
    if ~QUIET
        fprintf('%3d\t%10.4f\t%10.4f\n',...
            iter,TNN,gap)               ;
    end

end

x_opt  = reshape(x,size)          ;

end
