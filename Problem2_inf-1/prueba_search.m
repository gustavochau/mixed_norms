clc;
clear;
close all;

N = 1000;
M = 300;
% rng(3);
B = (rand(N,M)-0.5);
lambda =0.1;
shrink = @(u,ll) sign(u).*max(abs(u)-ll,0);


% % %% CVX solution to corroborate results
% cvx_begin
%     variable A_cvx(N,M)
%     expression row_norm(N)
%     for ii=1:N
%         row_norm(ii) = norm(A_cvx(ii,:),1);
%     end
%     Mixed_norm = max(row_norm);
%     minimize(0.5*sum_square( B(:) - A_cvx(:)) + lambda* Mixed_norm )
% cvx_end



% possible_t = linspace(tau_1*lambda,norma_B*lambda,50);
% [ A_prueba, costo, costo0,tau_opt ] = solve_l1_search( B,lambda,possible_t);
% plot(possible_t*lambda,costo)
% max(abs(A_prueba(:)-A_cvx(:)))

tic
norma_B = max(sum(abs(B),2));
% possible_t = (0:0.01:norma_B)*lambda;
% for jj=1:N
%     tau_1_l(jj) = norm(shrink(B(jj,:),lambda),1);
% end
% [tau_1,i_tau] = max(tau_1_l);
[~,ib]=max(sum(abs(B),2));
% i_tau
% ib
tau_1 = norm(shrink(B(ib,:),lambda),1);
[ A_stef, tau_opt_stef,iter] = solve_l1_search_steffenson( B,lambda, tau_1);
% max(abs(A_stef(:)-A_cvx(:)))
X_stef = B-A_stef;
toc
abs(compute_mixed_norm(X_stef,1,inf)-lambda)

tic
[ X_sra, theta_opt ] = solve_sra(B,lambda);
toc
abs(compute_mixed_norm(X_sra,1,inf)-lambda)

% n = 0:0.1:4;
% c=0*n;
% for ii=1:length(n)
%     c(ii) = search_function_l1(n(ii),B,lambda);
% end
% plot(n,c)