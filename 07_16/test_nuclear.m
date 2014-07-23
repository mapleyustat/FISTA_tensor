clear all
close all
clc

addpath('proxFunctions/')                        ;

%% problem set up

n1             =  20                             ;
n2             =  20                             ;
r              =  5                              ;
p              =  0.5                            ;

X_true         =  randn(n1,r)*randn(r,n2)        ;
x_true         =  X_true(:)                      ;
A              =  randn(round(p*n1*n2),n1*n2)    ;
b              =  A*x_true                       ;

lambda         =  1                              ;
L              =  norm(A)^2                      ;   % Lipschitz constant.
                                                     % equal to max(eig(A'*A))

%% function handles

f              =  @(x)1/2*norm(b-A*x)^2          ;
g              =  @(x)lambda*sum(svd(...
                           reshape(x,[n1,n2])))  ;
F              =  @(x)f+g                        ;
Grad_f         =  @(x)A'*(A*x-b)                 ; % Gradient operator of f
Prox_g         =  @(x,tau)proxF_nuc(x,...
                              lambda*tau,[n1,n2]); % Proximal function of g

%% FISTA
options.maxite =  1000                           ;
options.quiet  =  0                              ;
options.err    =  1e-2                           ;

x              =  randn(size(x_true))            ;
opt            =  FISTA_const(x,f,g,lambda,...
                  Prox_g,Grad_f,L,options)       ;
x_result       =  opt.x(:,options.maxite)        ;
X_result       =  reshape(x_result,[n1,n2])      ;

figure;
subplot(121) ; imshow(X_true)   ;title('x true') ;
subplot(122) ; imshow(X_result) ;title('result') ;