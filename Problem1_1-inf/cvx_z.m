function [ z ] = cvx_z( y,rho,x,G,v,num_groups,N ,lambda)
%CVX_X Summary of this function goes here
%   Detailed explanation goes here
        %% x update
    cvx_begin quiet
        variable z(num_groups,N)
        Mixed_norm = 0;
        for ii=1:num_groups
            Mixed_norm = Mixed_norm + norm(z(ii,:),Inf);
        end
        minimize( rho*0.5*sum_square( vec(z-G*diag(x)+v)) + lambda*Mixed_norm)
    cvx_end

end

