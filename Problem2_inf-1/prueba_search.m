clc;
clear;
close all;

N = 50;
M = 10;
B = (rand(N,M)-0.5);
lambda = 0.1;
shrink = @(u,ll) sign(u).*max(abs(u)-ll,0);

% %% CVX solution to corroborate results
cvx_begin
    variable A_cvx(N,M)
    expression row_norm(N)
    for ii=1:N
        row_norm(ii) = norm(A_cvx(ii,:),1);
    end
    Mixed_norm = max(row_norm);
    minimize(0.5*sum_square( B(:) - A_cvx(:)) + lambda* Mixed_norm )
cvx_end

norma_B = max(sum(abs(B),2));
% possible_t = (0:0.01:norma_B)*lambda;
for jj=1:N
    tau_1_l(jj) = norm(shrink(B(jj,:),lambda),1);
end
tau_1 = max(tau_1_l);



possible_t = linspace(tau_1*lambda,norma_B*lambda,50);
[ A_prueba, costo, costo0 ] = solve_l1_search( B,lambda,possible_t);
plot(possible_t/lambda,costo)
hold on
plot(possible_t/lambda,costo0)
max(abs(A_prueba(:)-A_cvx(:)))

% compute_mixed_norm(B-A_prueba,1,inf)-lambda
