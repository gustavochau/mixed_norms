function [x_admm,stats,u_iter] = admm_1_inf_overlap(y,groups,lambda)

%     N = 5;
%     M = 4;
%     B = rand(N,M);
    N = length(groups); %nmumber of groups, number of rows
    M = length(y); %elements in vector, number of columns
    
    G = zeros(N,M);
    for ii=1:N
        G(ii,groups{ii})=1;
    end
    
    G_aux = [];
    for ii=1:N
        G_aux = [G_aux; spdiags(G(ii,:)', 0, M,M)];
    end    
    
    opts.maxiter = 2; % máximo número de iteraciones
    opts.rho0 = 1; % rho inicial
    opts.tol = [1E-4 1E-4]; % [tolerancia absoluta  tolerancia_relativa]
    opts.parrho = [5 1.5 1.5];
    opts.rhoopt = 'fix'; % opción de rho

    % %%%%%%%% general ADMM
    params.y = y;
    params.lambda = lambda;
    params.mat_form = @(w) reshape(w,[N,M]);
    params.G_aux=G_aux;
    
    opts.x0 = rand(M,1);
    opts.z0 = vec(G*diag(opts.x0));
    opts.verbose = 1;

    solvers{1} = @xupdate_1inf_overlap;
    solvers{2} = @zupdate_1inf_overlap;
    matrices{1} = sparse(G); %A
    matrices{2} = -1; %B
    matrices{3} = zeros(N*M,1); %c
    matrices{4} = sparse(G'); %A^T
    func_costo = @(x,z) 0.5 * sum((x-y).^2) + lambda*(compute_mixed_norm(reshape(z,[N,M]),1,inf));
    [u_iter,stats,loops] = generalADMM(matrices, func_costo, solvers,  params, opts);
    x_admm = u_iter(loops).x;

end

