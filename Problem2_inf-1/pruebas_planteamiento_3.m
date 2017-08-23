clc;
clear;
close all;

N = 5;
M = 3;



num_test = 1;

% for pp=1:num_test
% rng(3);
B = (rand(N,M));
% B = ones(N,M);

% B= [1  1 1; 0.8 0.8 0.8; 0.8 0.8 0.8]
% B(1,:) = B(1,:)*2;
lambda = 0.2;
shrink = @(u,ll) sign(u).*max(abs(u)-ll,0);
% %% CVX solution to corroborate results
cvx_begin quiet
variable A_cvx(N,M)
expression row_norm(N)
for ii=1:N
    row_norm(ii) = norm(A_cvx(ii,:),1);
end
Mixed_norm = max(row_norm);
minimize( 0.5*sum_square( B(:) - A_cvx(:)) + lambda* Mixed_norm )
cvx_end
% 
% 
% for ii=1:N
%     A_temp = zeros(size(B));
%     A_temp(ii,:) = shrink(B(ii,:));
%     tau =norm(A_temp(ii,:),1);
%     for jj=1:N
%         if jj~=ii
%                         A_temp(jj,:) = projL1Mich(B(jj,:)', tau, 20);
% %             cvx_begin
% %             variable alpha_cvx(M,1)
% %             minimize( 0.5*sum_square( B(jj,:)' - alpha_cvx))
% %             subject to
% %             norm(alpha_cvx,1) <= tau
% %             cvx_end
% %             A_temp(jj,:) = alpha_cvx;
%             
%         end
%         
%     end
%     costo(ii) = tau;  %0.5*(norm(B-A_temp,'fro')^2) + compute_mixed_norm(A_temp,inf,1);
%     error(ii) = compute_mixed_norm(B-A_temp,1,inf)-lambda;
% end
% 
% [~,ind] = max(costo);
% A_test(ind,:) = shrink(B(ind,:));
% tau =norm(A_test(ind,:),1);
% for jj=1:N
%     if jj~=ind
%                     A_temp(jj,:) = projL1Mich(B(jj,:)', tau, 20);
%         cvx_begin
%         variable alpha_cvx(M,1)
%         minimize( 0.5*sum_square( B(jj,:)' - alpha_cvx))
%         subject to
%         norm(alpha_cvx,1) <= tau
%         cvx_end
%         A_test(jj,:) = alpha_cvx;
%     end
% 
% end
% 
% sum(abs(A_cvx),2)
% 
% abs(A_test-A_cvx)
% 
for ii=1:N
    A_sur(ii,:) = shrink(B(ii,:),lambda);
    costo(ii) = norm(A_sur,1);
end
[tau,ind] = max(costo);

for jj=1:N
    A_temp(jj,:) = projL1Mich(B(jj,:)', tau, 20);
end
A_temp(ind,:) = shrink(B(ind,:),lambda);
% end

% disp('errors')
% % % mean(error)
% mean(error2)
% mean(error3)
