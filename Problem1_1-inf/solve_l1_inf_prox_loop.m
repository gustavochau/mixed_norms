function [ A, loops ] = solve_l1_inf_prox_loop( B, lambda, nMaxIter )
% SOLVE_L1_INF_REG Solves the mixed 1, inf norm regularization problem 
%   argmin_A {0.5||A-B||^2_F + lambda ||A||_{1,inf}},
% by means of separable l1-ball projections over the rows of B. 
% The l1-ball projections are performed with the accelerated newton's method 
% proposed by Rodriguez (2017).
%   B: input matrix
%   lambda: regularization parameter
%   nMaxIter: maximum number of iterationf for l1-ball projections
%   A: found optimal solution
%   loops: number of loops per row for the l1-ball projection

[nrows,ncols] = size(B);

loops = zeros(nrows,1);
A = zeros(nrows,ncols);

for ii=1:nrows
    [x, ~, loops(ii)] = projL1AccNewton(B(ii,:)', lambda, nMaxIter);
    A(ii,:) = B(ii,:) - x';
end

end

