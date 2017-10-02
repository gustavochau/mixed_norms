clc;
clear ;
close all;

addpath('../../../Problem2_inf-1');

data_folder = '';

% load mris
load([data_folder 'mris.mat']);
% load cooccurence
load([data_folder 'co_mat.mat']);

X=b;
clear b;
Y = mris';
% Y = Y/max(abs(Y(:)));
clear mris;

X=X(:,1:2:end);
Y = Y(:,1:4:end);

K = size(Y,2); %num tasks = voxels
N = size(Y,1); %num words
M = size(X,2); % num features

%% normalize
for ii=1:size(X,1)
    X(ii,:) = X(ii,:)/norm(X(ii,:),2);
end

% load('W_nr.mat')

W0 = 0*(rand(M,K)-0.5);
mem=0;
alpha=0.01;
max_iter=100;
gamma= 0.01;
tau = 2.5;
costo = @(U) costo_func(U,X,Y);
% costo = @(U) compute_mixed_norm(U,1,inf)-tau;
grad_op = @(U) gradient_MTL( Y,X,U );

proj_op = @(U,tau0) proj_newton_pruned( U, tau, tau0);
% proj_op = @(U,tau0) proj_mixed_general( U, tau, @prox12, 1,2,tau0);
dual_norm = @(U) compute_mixed_norm(U,Inf,1);
[W, num_iter ,hist,tiempo,hist_tau,cambio,h_bound,nw] = proj_grad_desc( W0, proj_op, grad_op, 1E-6, max_iter,alpha, 0, dual_norm);
[W_mem, ~ ,~,tiempo_mem,hist_tau_mem] = proj_grad_desc( W0, proj_op, grad_op, 1E-6, max_iter,alpha, 1, dual_norm);


figure
plot(cambio)
xlabel('# iteration k')
ylabel('||W_{(k)}-W_{(k-1)}||_{\infty,1}')

figure
semilogy(hist_tau_mem)
hold on
semilogy(0.99*h_bound)
% semilogy(ng*alpha)
% semilogy(nw)
legend('tau','lower bound')
xlabel('# iteration k')

