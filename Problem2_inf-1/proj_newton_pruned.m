function [ X_newt,tau_opt,iter ] = proj_newton_pruned( B, lambda )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
        N = size(B,1);
        for ii=1:N
            A_sur(ii,:) = shrink(B(ii,:),lambda);
            costo(ii) = norm(A_sur(ii,:),1);
        end
        [tau_1] = max(costo);
        [ A_newt, tau_opt,iter] = solve_l1_search_newton_pruned( B,lambda, tau_1);
        X_newt = B-A_newt;
end

