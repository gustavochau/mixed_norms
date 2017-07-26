clc;
clear;
close all;

N = 50;
M = 50;
B = (rand(N,M));
lambda = 0.01;

% %% CVX solution to corroborate results
cvx_begin
    variable A_cvx(N,M)
    expression row_norm(N)
    for ii=1:N
        row_norm(ii) = norm(A_cvx(ii,:),1);
    end
    Mixed_norm = max(row_norm);
    minimize( norm( B - A_cvx, 'fro' ) + lambda* Mixed_norm )
cvx_end

[A_admm,stats] = admm_inf_1(B,lambda);

mean(abs(A_cvx(:)-A_admm(:)))
plot(stats.cost)
% %% Loop method
% tic
% [ A, loops ] = solve_l1_inf_prox_loop( B, lambda, 20 );
% % error_cvx_loop = norm(A-A_cvx,'fro')/N/M;
% toc
% % block method
% 
% tic
% [ A, loops ] = solve_l1_inf_block( B, lambda, 20 );
% % error_cvx_block = norm(A-A_cvx,'fro')/N/M;
% toc