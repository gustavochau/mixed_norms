function [ z ] = cvx_z(x,rho,L,S,v,N,lambda)
%CVX_X Summary of this function goes here
%   Detailed explanation goes here
        %% x update
    cvx_begin quiet
        variable z(N,1)
        minimize( rho*0.5*sum_square( L*S*x-z+v) + lambda*norm(z,inf))
    cvx_end

end

