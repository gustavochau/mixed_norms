clc;
clear;
close all;

load('prueba.mat')

N = length(x);


x_func = mixed_xupdate1(x,z,A,B,c,v,rho,params, AT);

beta = params.beta;
L = params.L;
A = L*(diag(sign(x)));

cvx_begin
    variable x_cvx(N)
    minimize( 0.5*sum_square( x_cvx -beta ) + 0.5*rho* sum_square(A*x_cvx-z+v) )
cvx_end

max(abs(x_func - x_cvx))