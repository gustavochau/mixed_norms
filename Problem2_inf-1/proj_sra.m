function [ X_newt,tau_opt,iter ] = proj_sra( B, lambda, tau_0 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
        if compute_mixed_norm(B,1,inf)<lambda
            X_newt = B;
            return;
        end
        if nargin <3
            tau_0 = 0;
        end
        N = size(B,1);
%         for ii=1:N
%             A_sur(ii,:) = shrink(B(ii,:),lambda);
%             costo(ii) = norm(A_sur(ii,:),1);
% % %             costo(ii) = norm(B(ii,:),1);
%         end
%         [tau_1] = lambda*max(costo);
%         if tau_1 <= 0.0001
%             tau_1=0.01;
%         end
%         tau_1=0.01;
        [ X_newt, tau_opt,iter] = solve_sra( B,lambda,tau_0);
        
end

