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

X=X(:,1:end);
Y = Y(:,1:2:end);

K = size(Y,2); %num tasks = voxels
N = size(Y,1); %num words
M = size(X,2); % num features

%% normalize
for ii=1:size(X,1)
    X(ii,:) = X(ii,:)/norm(X(ii,:),2);
end

load('W_nr.mat')

W0 = W_nr; %0*(rand(M,K)-0.5);
alpha=0.02;
max_iter=500;
gamma= 0.01;
tau = compute_mixed_norm(W0,1,inf)*gamma;
costo = @(U) costo_func(U,X,Y);
% costo = @(U) compute_mixed_norm(U,1,inf)-tau;
grad_op = @(U) gradient_MTL( Y,X,U );

proj_op = @(U,tau0) proj_newton_pruned( U, tau, tau0);
[W_newton, num_iter_newton,hist_newton,tiempos_newton,costo_newton] = proj_grad_desc(W0, proj_op, grad_op, 1E-3, max_iter,alpha);
% 
% proj_op = @(U,tau0) proj_sra( U, tau,tau0);
% [W_sra, num_iter_sra,hist_sra,tiempos_sra,costo_sra] = proj_grad_desc(W0, proj_op, grad_op, 1E-3, max_iter,alpha,costo);
% 
% clear X Y U
% save(['mri_' num2str(tau) '.mat'])
