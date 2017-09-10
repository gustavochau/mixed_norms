function [ G ] = gradient_MTL( Y,X,W )
%GRADIENT_MTL Summary of this function goes here
%   Detailed explanation goes here
    [N,M] = size(X);
    K = size(Y,2);
    G = zeros(M,K);
    for ii=1:K
        G(:,ii) = X'*(X*W(:,ii)-Y(:,ii));
    end
end

