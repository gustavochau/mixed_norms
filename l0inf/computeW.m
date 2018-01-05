function [ W ] = computeW( X, tol)
%COMPUTEW Summary of this function goes here
%   Detailed explanation goes here
if nargin <2
    tol = 1E-4;
end
W = 1./X;
W(abs(X)<tol) = 0;

end

