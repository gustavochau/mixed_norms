clc;
clear all;
close all;

N = 100;
M = 100;

B = (rand(N,M)-0.5);
lambda = 0.1;

% %% CVX solution to corroborate results
cvx_begin 
    variable A_cvx(N,M)
    expression row_norm(N)
    for ii=1:N
        row_norm(ii) = norm(A_cvx(ii,:),1);
    end
    Mixed_norm = max(row_norm);
    minimize( 0.5*sum_square( B(:) - A_cvx(:)) + lambda* Mixed_norm )
cvx_end

all(sign(A_cvx(:)) == sign(B(:)))
