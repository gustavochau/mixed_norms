clc;
clear;
close all;

N = 5;
M = 5;
B = (rand(N,M));
lambda = 01;

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

jj=1;
% possible_t = (0:0.01:norma_B)*lambda;
possible_t = linspace(0*norma_B*lambda,1.1*norma_B*lambda,50);
for t = possible_t
    cvx_begin
    variable A_prueba(N,M)
    minimize( 0.5*sum_square( B(:) - A_prueba(:)) + t )
    subject to
        (sum(abs(A_prueba),2)) <= t/lambda;
    cvx_end
%     [A_prueba] = batch_projL1AccNewton(B, t/lambda, 10);
    costo(jj) = 0.5*sum_square( B(:) - A_prueba(:)) + lambda*max(sum(abs(A_prueba),2));
    a_store{jj} = A_prueba;
    jj = jj+1;
end

[~,ind] = min(costo);
A_prueba = a_store{ind};
plot(costo)
max(abs(A_prueba(:)-A_cvx(:)))
