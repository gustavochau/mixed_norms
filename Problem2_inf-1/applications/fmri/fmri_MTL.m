clc;
clear ;
close all;

addpath('/home/gustavo/Documents/MATLAB/Mixed_norms/Problem2_inf-1');

data_folder = '/home/gustavo/Gdrive/Papers/Mixed norm/';

% load mris
load([data_folder 'mris.mat']);
% load cooccurence
load([data_folder 'co_mat.mat']);

tau = 0.1;
X=b;
clear b;
Y = mris';
Y = Y/max(abs(Y(:)));
clear mris;

% X=X(:,1:10000);
% Y = Y(:,1:5000);

K = size(Y,2); %num tasks = voxels
N = size(Y,1); %num words
M = size(X,2); % num features

%% normalize
for ii=1:size(X,1)
    X(ii,:) = X(ii,:)/norm(X(ii,:),2);
end

% X = X./repmat(max(abs(X)),[size(X,1) 1]);

grad_op = @(U) gradient_MTL( Y,X,U );
% proj_op = @(U) proj_cvx( U, tau );
proj_op = @(U) proj_newton_pruned( U, tau );
% proj_op = @(U) U;
% rng(5);
W0 = (rand(M,K)-0.5);
[W, num_iter,hist,tiempos] = proj_grad_desc(W0, proj_op, grad_op, 1E-3, 1000,0.01);

