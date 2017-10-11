function [ A ] = prox_cvx( B, theta, norm2 )
%PROX_CVX Summary of this function goes here
%   Detailed explanation goes here
theta
cvx_clear 
    [N,M] = size(B);
    A = zeros(N,M);
%     for ii=1:N
%         cvx_begin quiet
%             variable row_a(1,M)
%             minimize( 0.5*sum_square( B(ii,:) - row_a) + theta*norm(row_a,norm2));
%         cvx_end
%        
%         
%         A(ii,:) = row_a;
%     end
    cvx_begin quiet
    variable A(N,M)
    expression row_norm(N)
    for ii=1:N
        row_norm(ii) = norm(A(ii,:),norm2);
    end
    Mixed_norm = sum(row_norm);
    minimize( 0.5*sum_square( B(:) - A(:)) +theta*(Mixed_norm)) 
    cvx_end
end

