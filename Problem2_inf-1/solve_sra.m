function [ X, theta_opt,iter] = solve_sra(Y,gamma)
%SOLVE_SRA Summary of this function goes here
%   Detailed explanation goes here
    if (compute_mixed_norm(Y,1,inf) <= gamma)
        X = Y;
        theta_opt = 0;
        iter = 0;
        return
    else
        f = @(theta) search_function_l1(theta,Y,gamma);
        options = optimset('TolX',1E-10);
        [theta_opt,~,~,output] = fzero(f,0,options);
        X = Y - loop_projL1Mich(Y, theta_opt, 20);
        iter = output.iterations;
    end

end

