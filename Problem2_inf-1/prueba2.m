clc;
clear;
close all;

N = 5;
M = 5;
B = (rand(N,M));
lambda = 0.1;

% %% CVX solution to corroborate results
% cvx_begin
%     variable A_cvx(N,M)
%     expression row_norm(N)
%     for ii=1:N
%         row_norm(ii) = norm(A_cvx(ii,:),1);
%     end
%     Mixed_norm = max(row_norm);
%     minimize(0.5*sum_square( B(:) - A_cvx(:)) + lambda* Mixed_norm )
% cvx_end
% 
% norma_B = max(sum(abs(B),2));

cvx_begin
    variable temp(N,M)
    expression row_norm(N)
    for ii=1:N
        row_norm(ii) = norm(temp(ii,:),inf);
    end
    Mixed_norm = sum(row_norm);
    minimize(0.5*sum_square( B(:) - temp(:)))
    subject to
        Mixed_norm <= lambda;    
cvx_end

A_cvx2 = B-temp;

temp = 0*B;
gamma = 0;
alpha=0.01;
for ii=1:500
    % minimize
        temp = solve_l1_inf_block( B, gamma, 20 );
    % update dual
        gamma = gamma + alpha*(  sum(max(abs(temp),[],2)) - lambda );
        costo(ii) = 0.5*sum_square( B(:) - temp(:)) + gamma*(sum(max(abs(temp),[],2))-lambda);
        gap(ii) = norm(B-temp-A_cvx2,'fro');
        
        delta(ii)=sum(max(abs(temp),[],2))-lambda;
end

A_ball = B-temp;

max(abs(A_ball(:)-A_cvx2(:)))
