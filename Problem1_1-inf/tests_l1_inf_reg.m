clc;
clear;
close all;

N = 1000;
M = 100;
B = rand(N,M);
lambda = 0.1;

% %% CVX solution to corroborate results
% cvx_begin
%     variable A_cvx(N,M)
%     Mixed_norm = 0;
%     for ii=1:N
%         Mixed_norm = Mixed_norm + norm(A_cvx(ii,:),Inf);
%     end
%     minimize( norm( B - A_cvx, 'fro' ) + lambda* Mixed_norm )
% cvx_end

%% Loop method
tic
[ A, loops ] = solve_l1_inf_prox_loop( B, lambda, 20 );
% error_cvx_loop = norm(A-A_cvx,'fro')/N/M;
toc
% block method

tic
[ A, loops ] = solve_l1_inf_block( B, lambda, 20 );
% error_cvx_block = norm(A-A_cvx,'fro')/N/M;
toc