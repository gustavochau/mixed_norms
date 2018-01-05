function [ A_newt,tau_opt,iter ] = prox_newton_pruned( B, lambda, tau0 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
        if compute_mixed_norm(B,1,inf)<lambda
            X_newt = 0;
            return;
        end
        
        if nargin < 3 || tau0==0
        
            N = size(B,1);
            if all(abs(B(:))<lambda)
               tau_1=0.0001;
            else
               tau_1=0;
               for ii=1:N
                   if any(abs(B(ii,:))>lambda)
                            tt = norm(shrink(B(ii,:),lambda),1);
                            if tt> tau_1
                                tau_1 = tt;
                            end
                   end;
               end    
            end

        else
            tau_1 = tau0;
        end
        
        [ A_newt, tau_opt,iter] = solve_l1_search_newton_pruned( B,lambda, tau_1);
end

