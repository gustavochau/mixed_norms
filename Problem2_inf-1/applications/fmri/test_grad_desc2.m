clc;
clear ;
close all;

addpath('/home/gustavo/Documents/MATLAB/mixed_norms/Problem2_inf-1');

K = 20; %num tasks = voxels
N = 30; %num words
M = 100; % num features
rng(45)
Y =3*(randn(N,K));
X = (randn(N,M));

tau = 1.5;
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

proj_op = @(U,tau0) proj_newton_pruned( U, tau,tau0 );
W0 = 0*(rand(M,K)-0.5);
% [W, num_iter ] = proj_grad_desc(W0, proj_op, grad_op, 1E-4, 1000,0.1);
mem = 0;
max_iter = 1200;
[W, num_iter ,hist,tiempo,hist_tau,hist_x,hist_g, hist_alpha] = proj_grad_desc_exo( W0, proj_op, grad_op, 1E-4, max_iter, mem);
max(abs(W(:)-W_cvx(:)))
save('prueba_mtl')

for ii=1:num_iter
    nw(ii) = compute_mixed_norm(hist_x(:,:,ii),inf,1);
    ng(ii) = compute_mixed_norm(hist_g(:,:,ii),inf,1);
    if ii>2
%         ratio = compute_mixed_norm(hist_x(:,:,ii-2),inf,1)/compute_mixed_norm(hist_x(:,:,ii-1),inf,1);
%         no(ii) = compute_mixed_norm(hist_x(:,:,ii-1)-alpha*hist_g(:,:,ii),inf,1)-(1/ratio)*compute_mixed_norm(hist_x(:,:,ii-1),inf,1);
        cambio(ii) = compute_mixed_norm(hist_x(:,:,ii-1)-hist_x(:,:,ii-2),inf,1);
%         no(ii) = compute_mixed_norm(hist_x(:,:,ii-1)-alpha*hist_g(:,:,ii),inf,1)-compute_mixed_norm(hist_x(:,:,ii-1),inf,1)-cambio(ii);
        no(ii) = compute_mixed_norm(hist_x(:,:,ii-1)-alpha*hist_g(:,:,ii),inf,1)- compute_mixed_norm(hist_x(:,:,ii-2)-alpha*hist_g(:,:,ii-1),inf,1);
    end
end

tt = [0,hist_tau(1:(num_iter-1))] - [0,alpha*ng(1:(num_iter-1))];
tt2 = [0,hist_tau(1:(num_iter-1))] + [0,alpha*ng(1:(num_iter-1))];

figure
semilogy(hist_tau)
hold on
semilogy(0.99*no)
semilogy(ng*alpha)
% semilogy(nw)
legend('tau','0.95no','nw')
