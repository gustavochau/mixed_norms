function [ A_prueba, tau_opt,iter] = solve_l1_search_newton_pruned( B,lambda, tau_0)
%SOLVE_L1_SEARCH Summary of this function goes here
%   Detailed explanation goes here
    norm0 = @(u) sum(abs(u)>0.0001,2);
    max_iter=1000;
    tol = 1E-10;
    nb = sum(abs(B),2);
    f_newt = @(tau) search_function_l1_der_pruned(tau,B,lambda,nb);
    [tau_opt,x_hist,error,iter]=root_newt(f_newt,tau_0,tol,max_iter);

%     options = optimset('TolX',1E-10);
%     [tau_opt,~,~,output] = fzero(f,tau_0,options);
    
    A_prueba = loop_projL1Mich(B, tau_opt, 20);
%     iter = output.iterations;
%     jj=1;
%     [N,M] = size(B);
%     costo = zeros(length(possible_t),1);
%     costo0 = zeros(length(possible_t),1);
%     for t = possible_t
%         [A_prueba] = loop_projL1Mich(B, t/lambda, 20);
%         costo(jj) = 0.5*sum_square( B(:) - A_prueba(:)) + lambda*max(sum(abs(A_prueba),2));
%         costo0(jj) = 0.5*sum_square( B(:) - A_prueba(:)) + lambda*max(norm0(A_prueba));
%         a_store{jj} = A_prueba;
%         jj = jj+1;
%     end
% 
% [tau_opt,ind] = min(costo);
% A_prueba = a_store{ind};

end

