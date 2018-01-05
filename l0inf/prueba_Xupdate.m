clc;
clear;
close all;

N = 5;
M = 2;
B = rand(N,M);
X = rand(N,M);
Z = rand(N,M);
V = rand(N,M);
W = 1./X;
W(abs(X)<1E-6) = 0;

rho=1;
cvx_begin
    variable X_cvx(N,M)
    minimize  (0.5*sum_square(X_cvx(:)-B(:)) +0.5*rho* sum_square(vec(W.*X_cvx-Z+V)))
cvx_end

X_grad = (rho*W.*(Z-V)+B)./(1+rho*W.*W);