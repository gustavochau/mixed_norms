clc;
clear;
close all;

N = 100;
M = 20;

num_test = 10;

for pp=1:num_test
pp
B = (rand(N,M)-0.5);
lambda = 0.1;

% %% CVX solution to corroborate results
cvx_begin quiet
    variable A_cvx(N,M)
    expression row_norm(N)
    for ii=1:N
        row_norm(ii) = norm(A_cvx(ii,:),1);
    end
    Mixed_norm = max(row_norm);
    minimize( 0.5*sum_square( B(:) - A_cvx(:)) + lambda* Mixed_norm )
cvx_end


opts.maxiter = 300; % máximo número de iteraciones
opts.rho0 = 1; % rho inicial
opts.tol = [1E-6 1E-4]; % [tolerancia absoluta  tolerancia_relativa]
opts.parrho = [5 1.5 1.5];
opts.rhoopt = 'fix'; % opción de rho
opts.lambda = lambda;
opts.verbose = 0;

% tic
% [u_iter,stats,loops] = my_admm_inf_1(B,lambda, opts);
% times_admm(pp) = toc;
% 
% x_admm = u_iter(loops).x;
% A_admm = reshape(x_admm,[M,N])';
% 
% error(pp) = max(abs(A_cvx(:)-A_admm(:)));
% num_loops(pp) = loops;

tic
[u_iter,stats,loops] = my_admm_inf_1_v2(B,lambda, opts);
times_admm2(pp) = toc;

x_admm = u_iter(loops).x;
A_admm = reshape(x_admm,[M,N])';

error2(pp) = max(abs(A_cvx(:)-A_admm(:)));
num_loops2(pp) = loops;



tic
[u_iter,stats,loops] = my_admm_inf_1_v3(B,lambda, opts);
times_admm3(pp) = toc;

x_admm = u_iter(loops).x;
A_admm = reshape(x_admm,[M,N])';

error3(pp) = max(abs(A_cvx(:)-A_admm(:)));
num_loops3(pp) = loops;


% plot(stats.cost)
end

disp('errors')
% % mean(error)
mean(error2)
mean(error3)

disp('timess')
% mean(times_admm)
mean(times_admm2)
mean(times_admm3)
