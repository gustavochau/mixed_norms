clc;
clear;
close all;

N = 100;
B = rand(N,N);
lambda = 0.1;

% %% CVX solution to corroborate results
cvx_begin
    variable A_cvx(N,N)
    Mixed_norm = 0;
    for ii=1:N
        Mixed_norm = Mixed_norm + norm(A_cvx(ii,:),Inf);
    end
    minimize( norm( B - A_cvx, 'fro' ) + lambda* Mixed_norm )
cvx_end

%% Our method

[ A, loops ] = solve_l1_inf_prox_loop( B, lambda, 20 );
% error_cvx = norm(A-A_cvx,'fro')/N/N
