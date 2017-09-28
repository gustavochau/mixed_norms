function [ X_newt,tau_opt,iter ] = proj_steffensen( B, lambda, tau_0 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
        if compute_mixed_norm(B,1,inf)<lambda
            X_newt = B;
            return;
        end
        
        if nargin < 3
        
            N = size(B,1);
            if all(abs(B(:))<lambda)
               tau_1=0;
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
            tau_1 = tau_0;
        end
        
        [ A_newt, tau_opt,iter] = solve_l1_search_steffenson( B,lambda, tau_1);
        X_newt = B-A_newt;
end

