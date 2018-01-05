function [ g, dg] = search_function_l1_der_pruned(tau,B,lambda,absB)
%SEARCH_FUNCTION_L1 Summary of this function goes here
%   Detailed explanation goes here
%     if tau==0
%         g = compute_mixed_norm(B,1,inf)-lambda;
%         dg = -size(B,1)/(size(B,2));
%     else
        ind_util = logical(absB>=tau);
        B = B(ind_util,:);
        [A,l] = loop_projL1Mich(B, tau, 20);
    %     nb = sum(abs(B),2);
    %     s = logical(nb<tau/lambda);
    %     na = sum(abs(A),2);
    %     g = sum(na(~s)-tau/lambda) + sum(abs(na(s)-nb(s)));
        g = compute_mixed_norm(B-A,1,inf)-lambda;

        [nrows,ncols] = size(B);
        nb= sum(abs(B),2);
        dg=0;
        for ii=1:nrows
            if nb(ii) >= tau
                  z =  sign(B(ii,:)').*double(abs(B(ii,:))'>l(ii));
                  dg = dg-1/(z'*z);  
            end
        end
%     end
end

