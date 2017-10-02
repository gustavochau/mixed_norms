function [ X,theta_opt,iter ] = proj_mixed_general( B, lambda, prox_op, norm1,norm2,tau_0)
%UNTITLED Solves argmin_x 0.5 || x-y||^2_F s.t. ||x||_{norm1,norm2} <=
%lambda
%   tau_0: initial point for root search
        if compute_mixed_norm(B,norm1,norm2)<lambda
            X = B;
            theta_opt = 0;
            iter = 0;
            return;
        end
        if nargin <6
            tau_0 = 0;
        end
        g = @(theta) compute_mixed_norm(prox_op(B,theta),norm1,norm2)-lambda;
%         options = optimset('TolX',1E-10);
%         [theta_opt,~,~,output] = fzero(g,tau_0,options);
    tol_u = 1E-6;
    tol = 1E-5;
    max_iter=50;
        [theta_opt,~,~,iter]=steff_amat(g,tau_0,tol,max_iter,tol_u);

        X=prox_op(B,theta_opt);
%         iter = output.iterations;

end

