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
max_iter=50;
gamma= 0.01;
tau = 2.5;
costo = @(U) costo_func(U,X,Y);
% costo = @(U) compute_mixed_norm(U,1,inf)-tau;
grad_op = @(U) gradient_MTL( Y,X,U );

proj_op = @(U,tau0) proj_newton_pruned( U, tau, tau0);
[W, num_iter ,hist,tiempo,hist_tau] = proj_grad_desc( W0, proj_op, grad_op, 1E-4, max_iter,alpha, 0);
[W_mem, ~ ,~,tiempo_mem,hist_tau_mem,cambio,h_bound,nw] = proj_grad_desc( W0, proj_op, grad_op, 1E-4, max_iter,alpha, 1);


figure
plot(nw)

figure
semilogy(hist_tau)
hold on
semilogy(0.99*h_bound)
% semilogy(ng*alpha)
% semilogy(nw)
legend('tau','0.95no','nw')
