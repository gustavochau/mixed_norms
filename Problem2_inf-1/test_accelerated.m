clc
clear;
close all;


N = 500;
M = 200;
B = (rand(N,M)-0.5);
lambda =0.1;
norma_B = max(sum(abs(B),2));
    [~,ib]=max(sum(abs(B),2));
    tau_0 = norm(shrink(B(ib,:),lambda),1);
% tau_range = tau_1:0.01:tau_1*2;
% 
%     
% for ii=1:length(tau_range)
%     f(ii) = search_function_l1(tau_range(ii),B,lambda);
%     acc_f(ii) = search_function_l1_acc(tau_range(ii),B,lambda);
% end

% plot(tau_range,f)
% hold on
% plot(tau_range,acc_f)

f = @(tau) search_function_l1(tau,B,lambda);
f_acc = @(tau) search_function_l1_acc(tau,B,lambda);

max_iter=1000;
tol_u = 1E-10;
tol = 1E-10;
tic
[tau_opt,x_hist,error,iter]=steff_amat(f,tau_0,tol,max_iter,tol_u);
toc
tic
[tau_opt_acc,x_hist_acc,error_acc,iter_acc]=steff_amat(f_acc,tau_0,tol,max_iter,tol_u);
toc