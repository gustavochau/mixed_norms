clc;
clear;
close all;

N = 200;
M = 20;
B = rand(N,M);

lambda = 2;
maxiter = 15000;
opts.maxiter = maxiter; % máximo número de iteraciones
opts.rho0 =5; % rho inicial
opts.tol = [1E-4 1E-2]; % [tolerancia absoluta  tolerancia_relativa]
opts.parrho = [5 1.5 1.5];
opts.rhoopt = 'fix'; % opción de rho
opts.lambda = lambda;
opts.verbose = 1;

[u_iter,stats,loops] = my_admm_inf_0(B,lambda, opts);
X_admm = u_iter(loops).X;