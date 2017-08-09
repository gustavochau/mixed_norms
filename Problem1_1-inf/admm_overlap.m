function [u_iter,stats,loops] = admm_overlap(y, groups, handle_mixed_norm , solver_mixed_norm, opts)
%MYADMM Resuelve el problema de regularizacion
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
N = length(groups); % number of rows = number of groups
M = length(y); % number of columns = number of elements in target variable
mat_form = @(w) reshape(w,[M,N])'; % function for reshaping vectorized variable ROW-WISE

% history variables
costo = zeros(opts.maxiter,1);
costo_aumentada = zeros(opts.maxiter,1);
rho_hist = zeros(opts.maxiter,1);
tiempo = zeros(opts.maxiter,1);
primal_res = zeros(opts.maxiter,1);
dual_res = zeros(opts.maxiter,1);

% variables and parameters for admm
x = zeros(M,1); % first primal variable
z = zeros(N*M,1); % second primal variable
v = zeros(N*M,1); % dual variable
rho = opts.rho0;
lambda = opts.lambda;
loops=opts.maxiter;

G = zeros(N,M); % G matrix whose (i,j) element is 1 if j-th coordinate in y belongs to the i-th group
for ii=1:N
    G(ii,groups{ii})=1;
end

H = []; % H matrix that vectorizes the operation in G in a ROW-WISE fashion
for ii=1:N
    H = [H; spdiags(G(ii,:)', 0, M,M)];
end    

func_costo = @(x1) 0.5* sum((x1-y).^2) + lambda*handle_mixed_norm(G*diag(x1));
aug_func_costo = @(x1,z1,rho1) 0.5* sum((x1-y).^2) + lambda*handle_mixed_norm(mat_form(z1))+ 0.5*rho1*sum((H*x1-z1).^2);


for k=1:opts.maxiter
    
    disp(['Iteración ' num2str(k) ' de ADMM'])
    
    % ==================================
    % Calcular valor de la función costo
    % ==================================
    costo(k) = func_costo(x);
    costo_aumentada(k) =aug_func_costo(x,z,rho);

%     if (opts.verbose)
%         disp(['Total: ' num2str( costo_aumentada(k))])
%     end    
    t0 = tic; % inicio de contador de tiempo
    
    % =================
    % x update
    % =================
    x = xupdate_1inf_overlap( x,z,v,y,rho,H);
    
    if (opts.verbose)
        disp(['x actualizado: ' num2str(aug_func_costo(x,z,rho))])
    end
    
    % =================
    % z update
    % =================
    z_anterior = z;
    z = solver_mixed_norm( x,v,rho,lambda,G,mat_form);
    if (opts.verbose)
        disp(['z actualizado: ' num2str(aug_func_costo(x,z,rho))])
    end
    
    % ============================================
    % update of scaled dual variable
    % ============================================
    v = v + H*x-z;
    
    if (opts.verbose)

        disp(['v actualizado: ' num2str(aug_func_costo(x,z,rho))])
    end
    % ============================================
    % Stopping criteria
    % ============================================    
    
    primal_res(k) = norm(H*x-z,2);
    dual_res(k) = norm(rho*H'*(z-z_anterior),2);
   
    eps_primal(k)= sqrt(N*M)*opts.tol(1)+opts.tol(2)*max([norm(H*x,2),norm(z,2)]); % tolerancia para residuo primal
    eps_dual(k)=sqrt(M)*opts.tol(1)+opts.tol(2)*norm(H'*rho*v,2); % tolerancia para residuo dual
    
    tiempo(k) = toc(t0);
    u_iter(k).x = x;
    u_iter(k).z = z;
    
    if ((primal_res(k) <= eps_primal(k))  && (dual_res(k)<=eps_dual(k)))
        disp('Se alcanzaron condiciones de parada')
        loops = k;
        break;
    end
    
    % ========================================
    % Adaptar rho para la siguiente iteracción
    % ========================================
    
    if (strcmp(opts.rhoopt,'var'))
        mu = opts.parrho(1);
        tau_inc = opts.parrho(1);
        tau_dec = opts.parrho(2);
        
        if (primal_res> mu*dual_res)
            rho = tau_inc*rho; % aumentar rho
            v = v/tau_inc;  % reescalar v
        elseif (dual_res> mu*primal_res)
            rho = rho/tau_dec; % disminuir rho
            v = v*tau_dec; % reescalar v
        end
    end
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
stats.eps_primal=eps_primal;
stats.eps_dual = eps_dual;
end
