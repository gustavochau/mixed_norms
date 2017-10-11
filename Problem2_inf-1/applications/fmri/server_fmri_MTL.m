setenv('LC_ALL','C')

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

W0 = W_nr; % 0*(rand(M,K)-0.5); %W_nr; 
alpha=0.02;
max_iter=300;
gamma= 0.2;
tau = compute_mixed_norm(W0,1,inf)*gamma;
costo = @(U) costo_func(U,X,Y);
grad_op = @(U) gradient_MTL( Y,X,U );
dual_norm = @(U) compute_mixed_norm(U,Inf,1);

%proj_op = @(U,tau0) proj_newton_pruned( U, tau, tau0);
%[W_newton, num_iter_newton,hist_newton,tiempos_newton,cambio_newton,h_bound_newton] = proj_grad_desc( W0, %proj_op, grad_op, 1E-4, max_iter,alpha, 1, dual_norm);

%rows = max(abs(W_newton),[],2);
%save(['mri_' num2str(tau) '_newton_bound.mat'])
%save(['mri_' num2str(tau) '_newton_bound.mat'],'rows','num_iter_newton','hist_newton','tiempos_newton','cambio_newton','h_bound_newton');

%proj_op = @(U,tau0) proj_sra( U, tau,tau0);
%[W_sra, num_iter_sra,hist_sra,tiempos_sra,costo_sra,h_tau_sra] = proj_grad_desc( W0, proj_op, grad_op, 1E-4, max_iter,alpha, 1, dual_norm);

%rows = max(abs(W_sra),[],2);

%clear X Y U
%save(['mri_' num2str(tau) '_sra_bound.mat'],'rows','num_iter_sra','hist_sra','tiempos_sra','costo_sra','h_tau_sra');

proj_op = @(U,tau0) proj_steffensen_pruned( U, tau,tau0);
[W_ste, num_iter_ste,hist_ste,tiempos_ste,costo_ste,h_tau_ste] = proj_grad_desc( W0, proj_op, grad_op, 1E-4, max_iter,alpha, 1, dual_norm);

rows = max(abs(W_ste),[],2);

clear X Y U
save(['mri_' num2str(tau) '_ste_bound.mat'],'rows','num_iter_ste','hist_ste','tiempos_ste','costo_ste','h_tau_ste');
