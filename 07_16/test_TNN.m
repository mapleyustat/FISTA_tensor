clear all
close all
clc

addpath('proxFunctions/')                        ;
addpath('toolbox/')                              ;

%% problem set up

% n1             =  5                              ;
% n2             =  5                              ;
% n3             =  10                             ;
% Tsize          =  [n1,n2,n3]                     ;
% r              =  1                              ;
% p              =  0.5                            ;
% 
% X_true         =  randn(n1,n2,n3)                ;
% [U S V]        =  ntsvd(X_true)                  ;
% S(r+1:end,...
%     r+1:end,:) =  0                              ;
% X_true         =  tprod(tprod(U,S),tran(V))      ;
% x_true         =  X_true(:)                      ;

%% video data

load('basketball.mat') ;
video          =  imresize(vidd,0.04)             ;
nframes        =  10                             ;
X_true         =  video(:,:,1:nframes)           ;
[n1,n2,n3]     =  size(X_true)                   ;
Tsize          =  [n1,n2,n3]                     ;
p              =  0.3                            ;
r              =  3                              ;
[U S V]        =  ntsvd(X_true)                  ;
S(r+1:end,...
    r+1:end,:) =  0                              ;
X_true         =  tprod(tprod(U,S),tran(V))      ;
x_true         =  X_true(:)                      ;
%%

A              =  randn(round(p*n1*n2*n3),...
                                       n1*n2*n3) ;
b              =  A*x_true                       ;

lambda         =  1e6                           ;
L              =  norm(A)^2                      ;   % Lipschitz constant.
                                                     % equal to max(eig(A'*A))

%% function handles

f              =  @(x)1/2*norm(b-A*x)^2          ;
g              =  @(x)lambda*TNN(...
                           reshape(x,[n1,n2,n3]));
F              =  @(x)f+g                        ;
Grad_f         =  @(x)A'*(A*x-b)                 ; % Gradient operator of f
Prox_g         =  @(x,tau)proxF_TNN(x,...
                              lambda*tau,Tsize   ); % Proximal function of g

%% FISTA
options.maxite =  1000                           ;
options.quiet  =  0                              ;
options.err    =  1                               ;

x              =  randn(size(x_true))            ;
opt            =  FISTA_const(x,f,g,lambda,...
                  Prox_g,Grad_f,L,options)       ;
x_result       =  real(opt.x(:,options.maxite))  ;
X_result       =  reshape(x_result,Tsize)        ;

figure;
subplot(231) ; imagesc(X_true(:,:,1))   ;colormap(gray);title('x true1') ;
subplot(232) ; imagesc(X_true(:,:,2))   ;colormap(gray);title('x true2') ;
subplot(233) ; imagesc(X_true(:,:,3))   ;colormap(gray);title('x true3') ;
subplot(234) ; imagesc(X_result(:,:,1)) ;colormap(gray);title('result1') ;
subplot(235) ; imagesc(X_result(:,:,2)) ;colormap(gray);title('result2') ;
subplot(236) ; imagesc(X_result(:,:,3)) ;colormap(gray);title('result3') ;

figure;
for i= 1:n3
    subplot(121);imagesc(X_true(:,:,i));colormap(gray);title('true');
    subplot(122);imagesc(X_result(:,:,i));colormap(gray);title('result');
    pause(.1)
end
