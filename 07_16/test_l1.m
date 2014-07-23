clear all
close all
clc

addpath('proxFunctions/')                       ;

%% problem set up

n              =  600                           ;
p              =  0.3                           ;

A              =  randn(round(p*n),n)           ;
sign           =  rand(n,1)                     ;
x_true         =  randn(n,1).*(sign>0.95)       ;
b              =  A*x_true                      ;

lambda         =  1                             ;
L              =  norm(A)^2                     ;   % Lipschitz constant.
                                                   % equal to max(eig(A'*A))

%% function handles

f              =  @(x)1/2*norm(b-A*x)^2         ;
g              =  @(x)lambda*norm(x,1)          ;
F              =  @(x)f+g                       ;
Grad_f         =  @(x)A'*(A*x-b)                ; % Gradient operator of f
Prox_g         =  @(x,tau)proxF_l1(x,lambda*tau); % Proximal function of g

%% FISTA
options.maxite =  1000                          ;
options.quiet  =  0                             ;
options.err    =  1e-2                          ;

x              =  randn(size(x_true))           ;
opt            =  FISTA_const(x,f,g,lambda,...
                  Prox_g,Grad_f,L,options)      ;
x_result       =  opt.x(:,options.maxite)       ;

figure;
subplot(121); plot(x_true)   ; title('x true')  ;
subplot(122); plot(x_result) ; title('result')  ;