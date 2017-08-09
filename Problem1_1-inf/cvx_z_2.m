function [ z_v ] = cvx_z_2( rho,x,G,v,num_groups,N ,lambda)
%CVX_X Summary of this function goes here
%   Detailed explanation goes here
        %% x update
    cvx_begin quiet
        variable z(num_groups,N)
        Mixed_norm = 0;
        for ii=1:num_groups
            Mixed_norm = Mixed_norm + norm(z(ii,:),Inf);
        end
        minimize( rho*0.5*sum_square(vec(G*diag(x) -z + reshape(v,N,num_groups)')) + lambda*Mixed_norm)
    cvx_end

    z_v = vec(z');
end

