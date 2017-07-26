function [A_admm,stats] = admm_inf_1(B,lambda)

%     N = 5;
%     M = 4;
%     B = rand(N,M);
    [N,M] = size(B);
%     lambda = 0.01;

    opts.maxiter = 20; % máximo número de iteraciones
    opts.rho0 = 1; % rho inicial
    opts.tol = [1E-4 1E-4]; % [tolerancia absoluta  tolerancia_relativa]
    opts.parrho = [5 1.5 1.5];
    opts.rhoopt = 'var'; % opción de rho

    % %%%%%%%% general ADMM
    num_ele = N*M;
    B_t = B';
    params.beta = B_t(:);
    params.lambda = lambda;
    params.L = kron(speye(N),ones(1,M));
    opts.x0 = abs(rand(num_ele,1));
    opts.z0 = params.L*abs(opts.x0);
    opts.verbose = 1;

    solvers{1} = @mixed_xupdate1;
    solvers{2} = @mixed_zupdate1;
    matrices{1} = params.L*diag(sign(opts.x0)); %A
    matrices{2} = -speye(N); %B
    matrices{3} = zeros(N,1); %c
    matrices{4} = (params.L*sign(opts.x0))'; %A^T
    func_costo = @(x,z) 0.5 * sum((x-B_t(:)).^2) + lambda*(norm(z,inf));
    [u_iter,stats] = generalADMM(matrices, func_costo, solvers,  params, opts);
    A_admm = (reshape(u_iter(end).x,M,N))';

end

