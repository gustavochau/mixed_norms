function [ X, theta_opt ] = solve_sra(Y,gamma)
%SOLVE_SRA Summary of this function goes here
%   Detailed explanation goes here
    if (compute_mixed_norm(Y,1,inf) <= gamma)
        X = Y;
        theta = 0;
        iters = 0;
        return
    else
        f = @(theta) search_function_l1(theta,Y,gamma);
        options = optimset('TolX',1E-10);
        theta_opt = fzero(f,0,options);
        X = Y - loop_projL1Mich(Y, theta_opt, 20);
    end

end

