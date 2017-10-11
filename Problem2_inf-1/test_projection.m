clc;
clear;
close all;

N = 20;
M = 5;
B = (rand(N,M));
% theta = 0.2;
lambda = 3.5;
norm1 = 1;
norm2 = 2;
%% CVX solution to corroborate results
cvx_begin
    variable A_cvx(N,M)
    expression row_norm(N)
    for ii=1:N
        row_norm(ii) = norm(A_cvx(ii,:),norm2);
    end
    Mixed_norm = sum(row_norm);
    minimize( 0.5*sum_square( B(:) - A_cvx(:)) )
    subject to 
        Mixed_norm<= lambda;
cvx_end
% 
% cvx_begin
%     variable A_cvx(N,M)
%     expression row_norm(N)
%     for ii=1:N
%         row_norm(ii) = norm(A_cvx(ii,:),norm2);
%     end
%     Mixed_norm = sum(row_norm);
%     minimize( 0.5*sum_square( B(:) - A_cvx(:)) + theta* Mixed_norm )
% cvx_end

% prox = @(B,theta) prox_cvx( B, theta, norm2 );
[ A,theta_opt,iter ] = proj_mixed_general( B, lambda, @prox12, norm1,norm2,0);
% [ A,theta_opt,iter ] = proj_sra( B, lambda, 0 );
A = prox12(B,theta_opt);
max(abs(A(:)-A_cvx(:)))

compute_mixed_norm(A,norm1,norm2)
compute_mixed_norm(B-A,1/(1-1/norm1),1/(1-1/norm2))


