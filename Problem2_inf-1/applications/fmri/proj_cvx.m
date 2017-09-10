function [ A ] = proj_cvx( B, tau )
%PROJ_CVX Summary of this function goes here
%   Detailed explanation goes here
[M,K] = size(B);
cvx_begin 
    variable A(M,K)
    expression row_norm(M)
    for ii=1:M
        row_norm(ii) = norm(A(ii,:),inf);
    end
    Mixed_norm = sum(row_norm);
    minimize( 0.5*sum_square( B(:) - A(:) ))
    subject to
        Mixed_norm<=tau
cvx_end

end

