clc;
clear;
close all;

N = 5;
M = 3;
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
norma_B = max(sum(abs(B),2));
% possible_t = (0:0.01:norma_B)*lambda;
for jj=1:N
    tau_1_l(jj) = norm(shrink(B(jj,:),lambda),1);
end
tau_1 = max(tau_1_l);
possible_t_l1 = linspace(tau_1*lambda,norma_B*lambda,50);


%%% norm 0 limits
norma_B = max(norm0(B));
for jj=1:N
    tau_1_l(jj) = norm0(hard_thresh(B(jj,:),sqrt(2*lambda)));
end
tau_1 = max(tau_1_l);
possible_t_l2 = linspace(tau_1*lambda,1.3*norma_B*lambda,50);

possible_t = sort(union(possible_t_l2,possible_t_l1));

[ ~, costo_l1, costo_l1_0 , t_1] = solve_l1_search( B,lambda,possible_t);

[ A_prueba, costo_l0_pure, t_0] = solve_l0_search( B,lambda,possible_t);
plot(possible_t/lambda,costo_l0_pure)
hold on
% plot(possible_t/lambda,costo_l1)
plot(possible_t/lambda,costo_l1_0)
legend('l0-pure','solved w/l1 but cost in l0')
% max(abs(A_prueba(:)-A_cvx(:)))
t_0/t_1
% compute_mixed_norm(B-A_prueba,1,inf)-lambda
