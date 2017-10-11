function [ A ] = prox12( B, theta )
%PROX12 Summary of this function goes here
%   Detailed explanation goes here
    shrink = @(u,l) sign(u).*max(abs(u)-l,0);
    [N,M] = size(B);
    A = zeros(N,M);
    for ii=1:N
        nb = norm(B(ii,:),2);
        A(ii,:) = shrink(nb,theta)*B(ii,:)/nb;
    end
end

