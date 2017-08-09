function [x_admm,stats,u_iter] = admm_1_inf_overlap(y,groups,lambda,maxiter)

   
    opts.maxiter = maxiter; % máximo número de iteraciones
    opts.rho0 = 1; % rho inicial
    opts.tol = [1E-4 1E-2]; % [tolerancia absoluta  tolerancia_relativa]
    opts.parrho = [5 1.5 1.5];
    opts.rhoopt = 'fix'; % opción de rho
    opts.lambda = lambda;
    opts.verbose = 1;

    handle_mixed_norm = @(W) compute_mixed_norm(W,1,inf);
    solver_mixed_norm = @zupdate_1inf_overlap;
    
    [u_iter,stats,loops] = admm_overlap(y, groups, handle_mixed_norm , solver_mixed_norm, opts);
    x_admm = u_iter(loops).x;

end

