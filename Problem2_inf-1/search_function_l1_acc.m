function [ g ] = search_function_l1_acc(tau,B,lambda)
%SEARCH_FUNCTION_L1 Summary of this function goes here
%   Detailed explanation goes here
    A = loop_projL1Mich(B, tau, 20);
%     nb = sum(abs(B),2);
%     s = logical(nb<tau/lambda);
%     na = sum(abs(A),2);
%     g = sum(na(~s)-tau/lambda) + sum(abs(na(s)-nb(s)));
    g = tau*(compute_mixed_norm(B-A,1,inf)-lambda);

end

