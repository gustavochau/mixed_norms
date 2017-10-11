function [ A_prueba, tau_opt,iter] = solve_l1_search_steffenson_pruned( B,lambda, tau_0)
%SOLVE_L1_SEARCH Summary of this function goes here
%   Detailed explanation goes here
    norm0 = @(u) sum(abs(u)>0.0001,2);
    max_iter=1000;
    tol_u = 1E-6;
    tol = 1E-10;
    nb = sum(abs(B),2);
    f = @(tau) search_function_l1_pruned(tau,B,lambda,nb);
    [tau_opt,x_hist,error,iter]=steff_amat(f,tau_0,tol,max_iter,tol_u);

%     options = optimset('TolX',1E-10);
%     [tau_opt,~,~,output] = fzero(f,tau_0,options);
    
%     A_prueba = loop_projL1Mich(B, tau_opt, 20);
    A_prueba = loop_projL1AccNewt(B, tau_opt, 20);


end

