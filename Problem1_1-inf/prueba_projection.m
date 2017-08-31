clc;
clear;
close all;

N = 1000;
% rng(15)

b = rand(N,1)-0.5;
tau = 0.1;

% rng
cvx_begin 
    variable a_cvx(N,1)
    minimize( 0.5*sum_square( b - a_cvx))
    subject to
        norm(a_cvx,1) <= tau
cvx_end


[a_mich,lam_mich,l_mich] = projL1Mich(b, tau, 20);
[a_newt,lam_newt,l_newt] = projL1AccNewton_rfixed(b, tau, 20);

max(abs(a_mich - a_cvx))

max(abs(a_newt - a_cvx))