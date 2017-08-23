clc;
clear;
close all;

N = 30;
M = 20;
B = (rand(N,M)-0.5);
lambda = 0.1;
shrink = @(u,ll) sign(u).*max(abs(u)-ll,0);
hard_thresh = @(u,ll) u.*(abs(u)>ll);
norm0 = @(u) sum(abs(u)>0.0001,2);

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

%%% norm 1 limits


%%% norm 0 limits
norma_B = max(norm0(B));
for jj=1:N
    tau_1_l(jj) = norm0(hard_thresh(B(jj,:),sqrt(2*lambda)));
end
tau_1 = max(tau_1_l);

possible_t = linspace(tau_1*lambda,1.5*norma_B*lambda,50);
[ A_prueba, costo ] = solve_l0_search( B,lambda,possible_t);
plot(possible_t/lambda,costo)
% max(abs(A_prueba(:)-A_cvx(:)))

% compute_mixed_norm(B-A_prueba,1,inf)-lambda
