function [A_admm,stats] = admm_inf_0(B,lambda,maxiter)

    opts.maxiter = maxiter; % máximo número de iteraciones
    opts.rho0 = 1; % rho inicial
    opts.tol = [1E-4 1E-2]; % [tolerancia absoluta  tolerancia_relativa]
    opts.parrho = [5 1.5 1.5];
    opts.rhoopt = 'fix'; % opción de rho
    opts.lambda = lambda;
    opts.verbose = 1;

    [u_iter,stats,loops] = my_admm_inf_1(B,lambda, opts);
    x_admm = u_iter(loops).x;
end

