clc;
clear ;
close all;

addpath('/home/gustavo/Documents/MATLAB/mixed_norms/Problem2_inf-1');

K = 20; %num tasks = voxels
N = 60; %num words
M = 40; % num features
rng(45)
Y = 3*(rand(N,K) -0.5);
X = (rand(N,M) -0.5);

tau = 0.2;
cvx_begin 
    variable W_cvx(M,K)
    expression row_norm(N)
    expression obj_va
    for ii=1:M
        row_norm(ii) = norm(W_cvx(ii,:),inf);
    end
    Mixed_norm = sum(row_norm);
    obj_va = 0;
    for kk=1:K
        obj_va = obj_va + 0.5*sum_square( Y(:,kk) - X*W_cvx(:,kk));
    end
    minimize( obj_va )
    subject to
        Mixed_norm<=tau
cvx_end

grad_op = @(U) gradient_MTL( Y,X,U );
% proj_op = @(U) proj_cvx( U, tau );

proj_op = @(U) proj_newton_pruned( U, tau );
W0 = 0*(rand(M,K)-0.5);
[W, num_iter ] = proj_grad_desc(W0, proj_op, grad_op, 1E-4, 1000,0.1);

max(abs(W(:)-W_cvx(:)))
