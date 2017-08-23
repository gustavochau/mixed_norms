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

jj=1;

possible_t = linspace(tau_1*lambda,norma_B*lambda,50);
for t = possible_t
%     cvx_begin
%     variable A_prueba(N,M)
%     minimize( 0.5*sum_square( B(:) - A_prueba(:)) + t )
%     subject to
%         (sum(abs(A_prueba),2)) <= t/lambda;
%     cvx_end
    [A_prueba] = loop_projL1Mich(B, t/lambda, 20);
    costo(jj) = 0.5*sum_square( B(:) - A_prueba(:)) + lambda*max(sum(abs(A_prueba),2));
    a_store{jj} = A_prueba;
    jj = jj+1;
end



[tau_opt,ind] = min(costo);
A_prueba = a_store{ind};
plot(costo)
max(abs(A_prueba(:)-A_cvx(:)))

compute_mixed_norm(B-A_prueba,1,inf)-lambda
