function [u_iter,stats,loops] = my_admm_inf_0(B,lambda, opts)
%MYADMM Solves 0.5*||A-B||_F^2 + lambda * ||A||_{inf,1}
% PARÁMETROS DE ENTRADA
% =====================
%      matrices: Celda con las matrices A,B, C o los handles a las funciones que implementan Ax, Bz
% Cost_function: handle que calcula la función costo     
%       solvers: Celda con los solvers para los subproblemas de x y z
%          opts: Opciones como:
%           opts.maxiter: máximo número de operaciones
%            opts.rhoopt: 'fix' para rho fijo y 'var' para rho variable
%              opts.rho0: rho inicial
%            opts.parrho: parámetros para rho variable en la forma [mu, tau_incr, tau_dec]
%               opts.tol: vector de toleracion de la forma [eps_absoluto eps_relativo]
%                opts.x0: solución inicial
%                opts.z0: solución inicial
%                opt.verbose
% PARÁMETROS DE SALIDA
% =====================
%    u_iter: x,z de cada iteracion
%     stats: Información de la solución. Contiene:
%             stats.cost: valor de la función costo en cada iteración
%             stats.time: tiempo que toma cada iteración
%
% ====================
% Inicializar z y v y otras variable
% ====================

%inicialization
[N,M] = size(B); % N number of groups
beta = B(:)
mat_form = @(w) reshape(w,[M,N])'; % function for reshaping vectorized variable ROW-WISE

% history variables
costo = zeros(opts.maxiter,1);
costo_aumentada = zeros(opts.maxiter,1);
rho_hist = zeros(opts.maxiter,1);
tiempo = zeros(opts.maxiter,1);
primal_res = zeros(opts.maxiter,1);
dual_res = zeros(opts.maxiter,1);

% variables and parameters for admm
X = rand(N,M); % first primal variable
Z = rand(N,M); % second primal variable
V = zeros(N,M); % dual variable
W = computeW(X);

rho = opts.rho0;
loops=opts.maxiter;


func_costo = @(x1,w1) 0.5* sum((x1(:)-beta).^2) + lambda*compute_mixed_norm(w1.*x1,inf,1);
aug_func_costo = @(x1,z1,w1,rho1) 0.5* sum((x1(:)-beta).^2) + lambda*compute_mixed_norm(Z,1,inf)+ 0.5*rho1*sum(vec(w1.*x1-Z).^2);


for k=1:opts.maxiter
    
    if (opts.verbose)
        disp(['Iteración ' num2str(k) ' de ADMM'])
    end
    
    % ==================================
    % Calcular valor de la función costo
    % ==================================
    costo(k) = func_costo(X,W);
    costo_aumentada(k) =aug_func_costo(X,Z,W,rho);

%     if (opts.verbose)
%         disp(['Total: ' num2str( costo_aumentada(k))])
%     end    
    t0 = tic; % inicio de contador de tiempo
    
    % =================
    % x update
    % =================
%     alpha=cvx_x(beta,rho,z,L,S,v,N,M);
    X = xupdate_0inf( X,Z,V,B,rho,W);
    if (opts.verbose)
        disp(['x actualizado: ' num2str(func_costo(X,W))])
    end
%     S = sparse(diag(sign(alpha)));
    W = computeW(X);

    % =================
    % z update
    % =================
    Z_anterior = Z;
%     z = cvx_z(alpha,rho,L,S,v,N,lambda);
    Z = zupdate_0inf( X,V,rho,lambda,W);
%     z = solver_mixed_norm( alpha,v,rho,lambda,G,mat_form);
    if (opts.verbose)
        disp(['z actualizado: ' num2str(func_costo(X,W))])
    end
    
    % ============================================
    % update of scaled dual variable
    % ============================================
    V = V + W.*X-Z;
    
    if (opts.verbose)

        disp(['v actualizado: ' num2str(func_costo(X,W))])
    end
    % ============================================
    % Stopping criteria
    % ============================================    
    
    primal_res(k) = norm(W.*X-Z,'fro');
    dual_res(k) = norm(rho*W.*(Z-Z_anterior),'fro');
%    
%     eps_primal(k)= sqrt(N*M)*opts.tol(1)+opts.tol(2)*max([norm(L*(abs(X)),2),norm(Z,2)]); % tolerancia para residuo primal
%     eps_dual(k)=sqrt(M)*opts.tol(1)+opts.tol(2)*norm(rho*S'*L'*V,2); % tolerancia para residuo dual
%     
    tiempo(k) = toc(t0);
    u_iter(k).X = X;
    u_iter(k).Z = Z;
%     
%     if ((primal_res(k) <= eps_primal(k))  && (dual_res(k)<=eps_dual(k)))
%         disp('Se alcanzaron condiciones de parada')
%         loops = k;
%         break;
%     end
%     
%     % ========================================
%     % Adaptar rho para la siguiente iteracción
%     % ========================================
%     
%     if (strcmp(opts.rhoopt,'var'))
%         mu = opts.parrho(1);
%         tau_inc = opts.parrho(1);
%         tau_dec = opts.parrho(2);
%         
%         if (primal_res> mu*dual_res)
%             rho = tau_inc*rho; % aumentar rho
%             V = V/tau_inc;  % reescalar v
%         elseif (dual_res> mu*primal_res)
%             rho = rho/tau_dec; % disminuir rho
%             V = V*tau_dec; % reescalar v
%         end
%     end
    rho_hist(k) = rho;
    
    
end

% ==================================
% Dar formato a variables de salida
% ==================================

stats.cost = costo;
stats.costo_aumentada = costo_aumentada;
stats.time = tiempo;
stats.rho = rho_hist;
stats.primal_res=primal_res;
stats.dual_res=dual_res;
% stats.eps_primal=eps_primal;
% stats.eps_dual = eps_dual;
stats.loops = loops;
end
