clc;
clear;
close all;

load('prueba.mat')

M = length(z);

x = rand(length(x),1);

z_func = mixed_zupdate1( x,z,A,B,c,v,rho,params, AT );

L = params.L;
A = L*(diag(sign(x)));
lambda = params.lambda;

cvx_begin
    variable z_cvx(M)
    minimize( 0.5*rho*sum_square( A*x - z_cvx ) + lambda* norm(z_cvx,inf) )
cvx_end

max(abs(z_func - z_cvx))