clear all
close all
clc

n1 = 5;
n2 = 5;
n3 = 4;
p  = 0.3;

m = round(p*n1*n2*n3) ;
x = randn(n1,n2,n3)   ;
A = randn(m,100) ;
b = A*x(:)            ;
QUIET = 0;

x_opt = FISTA_const( A , b , [n1 n2 n3] , QUIET);