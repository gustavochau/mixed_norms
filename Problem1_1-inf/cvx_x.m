function [ x ] = cvx_x( y,rho,z,G,v,N )
%CVX_X Summary of this function goes here
%   Detailed explanation goes here
        %% x update
    cvx_begin quiet
        variable x(N,1)
        minimize( 0.5*sum_square( y-x) + 0.5*rho*sum_square(vec(z-G*diag(x)+v)))
    cvx_end

end

