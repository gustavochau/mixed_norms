function [u_iter,stats] = myADMM(Bases, In, lambda, opts)
%MYADMM Resuelve el problema de Basis Pursuit Denoising mediante ADMM
% PARÁMETROS DE ENTRADA
% =====================
%     Bases: Matriz con los handles a las funciones que implementan el
%            operador directo y traspuesto
%        In: Imagen ruidosa de entrada (vectorizada)
%    lambda: Parámetro de regularización
%      opts: Opciones como:
%           opts.maxiter: máximo número de operaciones
%            opts.rhoopt: 'fix' para rho fijo y 'var' para rho variable
%              opts.rho0: rho inicial
%            opts.parrho: parámetros para rho variable en la forma [mu, tau_incr, tau_dec]
%               opts.tol: vector de toleracion de la forma [eps_absoluto eps_relativo]
%                opts.u0: solución inicial
%                
% PARÁMETROS DE SALIDA
% =====================
%    u_iter: solución sparse en cada iteracion
%     stats: Información de la solución. Contiene:
%             stats.cost: valor de la función costo en cada iteración
%             stats.time: tiempo que toma cada iteración

% ====================
% Inicializar z y v y otras variable
% ====================

x0 = opts.u0; %inicializar

u_iter = zeros(size(x0,1),opts.maxiter); % array donde se guarda
xi = zeros(size(In));
z0 = x0;
v0 = zeros(size(x0,1),1);
costo = zeros(opts.maxiter,1);
costo_aumentada = zeros(opts.maxiter,1);
rho_hist = zeros(opts.maxiter,1);
tiempo = zeros(opts.maxiter,1);

Phi = Bases{1}; % operador directo
Phi_T = Bases{2}; %operador traspuesto

x = x0; % primera variable primal
z = z0; % segunda variable primal
v = v0; % variable dual
rho = opts.rho0;


func_costo = @(u) 0.5 * (norm(Phi(u)-In,2))^2 + lambda*(norm(u,1));
shrinkage  = @(vector,thres) sign(vector).*max(abs(vector) - thres,0);
func_costo_aumentada = @(x1,z1,v1,rho1) 0.5 * (norm(Phi(x1)-In,2))^2 + lambda*(norm(z1,1))+0.5*rho1*(norm(x1-z1+v1,2))^2-(rho1/2)*norm(v1);


for k=1:opts.maxiter
    
    disp(['Iteración ' num2str(k) ' de ADMM'])
    
    % ==================================
    % Calcular valor de la función costo
    % ==================================
    costo(k) = func_costo(x);
    costo_aumentada(k) =func_costo_aumentada(x,z,v,rho);
    disp(['Total: ' num2str( costo(k))])
    
    t0 = tic; % inicio de contador de tiempo
    
    % =================
    % minimización en x
    % =================
    % Consiste en la resolución de un sistema lineal
    % Como la matriz es simétrica e ímplicita se usará cgs
    
        
    
    
    disp(['Antes de todo: ' num2str(func_costo_aumentada(x,z,v,rho))])
    lhs = @(x) Phi_T(Phi(x)) + rho*x; % handle al lhs de la ecuación
        b = Phi_T(In)-rho*(v-z); % Lado derecho de la ecuación
        x = cgs(lhs,b,1E-5,30,[],[],x);
%     lhs = @(x) (1/rho)*Phi(Phi_T(x)) + x; % handle al lhs de la ecuación
%     w = (Phi_T(In)-rho*(v-z));
%     xi = cgs(lhs,Phi(w),1E-5,30,[],[],xi);
%     x = w/rho - Phi_T(xi)/(rho^2);
    
    disp(['x actualizado: ' num2str(func_costo_aumentada(x,z,v,rho))])
    
    % =================
    % minimización en z
    % =================
    % shrinkage
    z_anterior = z;
    z = shrinkage(x+v, lambda/rho);
    
    disp(['z actualizado: ' num2str(func_costo_aumentada(x,z,v,rho))])

    
    % ============================================
    % Actualización de la variable dual escalada
    % ============================================
    
    v = v + x - z;
    
    disp(['v actualizado: ' num2str(func_costo_aumentada(x,z,v,rho))])

    % ============================================
    % Condiciones de parada
    % ============================================    
    
    primal_res = norm(x - z);
    dual_res = norm(rho*(z-z_anterior));
    
    eps_primal= sqrt(size(x,1))*opts.tol(1)+opts.tol(2)*max([norm(x,2) norm(z,2)]); % tolerancia para residuo primal
    eps_dual=sqrt(size(x,1))*opts.tol(1)+opts.tol(2)*norm(rho*v); % tolerancia para residuo dual
    
    tiempo(k) = toc(t0);
    u_iter(:,k)= x;
        z_iter(:,k)= z;

    if ((primal_res <= eps_primal)  && (dual_res<=eps_dual))
        disp('Se alcanzaron condiciones de parada')
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
stats.z_iter = z_iter;
end
