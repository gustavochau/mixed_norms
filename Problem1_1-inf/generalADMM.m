function [u_iter,stats,loops] = generalADMM(matrices, func_costo, solvers,  params, opts)
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

%inicializacion

A = matrices{1};
B = matrices{2};
c = matrices{3};
AT = matrices{4}; %operador traspuesto


if(~isa(A,'function_handle'))
    A = @(w) matrices{1}*w;
end

if(~isa(AT,'function_handle'))
    AT = @(w) matrices{4}*w;
end

if(~isa(B,'function_handle'))
    B = @(w) matrices{2}*w;
end

solveX = solvers{1};
solveZ = solvers{2};
mat_form=params.mat_form;
% u_iter = zeros(size(x0,1),opts.maxiter); % array donde se guarda

costo = zeros(opts.maxiter,1);
costo_aumentada = zeros(opts.maxiter,1);
rho_hist = zeros(opts.maxiter,1);
tiempo = zeros(opts.maxiter,1);

x = opts.x0; % primera variable primal
z = opts.z0; % segunda variable primal
v = zeros(size(c)); % variable dual
rho = opts.rho0;

func_costo_aumentada = @(x1,z1,v1,rho1) func_costo(x1,z1) + (rho1/2)*norm(vec(A(diag(x))')+B(z1)-c+v1,2) - (rho1/2)*norm(v1);
loops=opts.maxiter;

for k=1:opts.maxiter
    
    disp(['Iteración ' num2str(k) ' de ADMM'])
    
    % ==================================
    % Calcular valor de la función costo
    % ==================================
    costo(k) = func_costo(x,z);
    
    costo_aumentada(k) =func_costo_aumentada(x,z,v,rho);
        if (opts.verbose)
            disp(['Total: ' num2str( costo_aumentada(k))])
        end    
    t0 = tic; % inicio de contador de tiempo
    
    % =================
    % x update
    % =================
%     x_anterior=x;
    x = solveX(x,z,A,B,c,v,rho,params, AT);
%         norm(x-x_anterior)

    if (opts.verbose)
    disp(['x actualizado: ' num2str(func_costo_aumentada(x,z,v,rho))])
    end
    % =================
    % z update
    % =================
    z_anterior = z;
    z = solveZ(x,z,A,B,c,v,rho,params, AT);
    if (opts.verbose)
    disp(['z actualizado: ' num2str(func_costo_aumentada(x,z,v,rho))])
    end
    % ============================================
    % update of scales dual variables
    % ============================================
%         v_anterior = v;
%     A = @(w) params.L*(sign(x).*w);
    v = v + vec((A(diag(x))')) + B(z) - c;
%                 norm(v-v_anterior)
    if (opts.verbose)

        disp(['v actualizado: ' num2str(func_costo_aumentada(x,z,v,rho))])
    end
    % ============================================
    % Condiciones de parada
    % ============================================    
    
    primal_res = norm(vec(A(diag(x))') + B(z) - c);
    dual_res = norm(rho*AT(params.mat_form(B(z-z_anterior))));
    
    eps_primal= sqrt(numel(c))*opts.tol(1)+opts.tol(2)*max([norm(vec(A(diag(x))'),'fro') norm(B(z),'fro') norm(c,'fro')]); % tolerancia para residuo primal
    eps_dual=sqrt(numel(x))*opts.tol(1)+opts.tol(2)*norm(AT(mat_form(rho*v)),'fro'); % tolerancia para residuo dual
    
    tiempo(k) = toc(t0);
    u_iter(k).x = x;
    u_iter(k).z = z;
    

    
    if ((primal_res <= eps_primal)  && (dual_res<=eps_dual))
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
end
