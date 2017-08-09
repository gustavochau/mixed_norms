function [ x ] = cvx_x_2( y,rho,z,H,v,N )
%CVX_X Summary of this function goes here
%   Detailed explanation goes here
        %% x update
    cvx_begin quiet
        variable x(N,1)
        minimize( 0.5*sum_square(x-y) + 0.5*rho*sum_square(H*x-z+v))
    cvx_end

end

