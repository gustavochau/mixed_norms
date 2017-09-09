clc;
clear;
close all;




N = 100;
M = 50;
% rng(3);

shrink = @(u,ll) sign(u).*max(abs(u)-ll,0);

num_real = 20;


pp=1;
possible_lambda = 0.01:0.01:0.3;
for ll=1:length(possible_lambda)
     disp([num2str(ll) '/' num2str(length(possible_lambda))])

    lambda = possible_lambda(ll);
    
    for zz=1:num_real
        %     rng(5*zz)
        
        B = ((rand(N,M)-0.5));
        
        nb,index = sort(sum(abs(B),2),'descending');
        
%         for ii=1:N
%             A_sur(ii,:) = shrink(B(ii,:),lambda);
%             costo(ii) = norm(A_sur(ii,:),1);
%         end
%         nb = sum(abs(B),2);
%         [tau_1] = max(costo);
%         num_may = sum(nb>tau_1);
%         for ii=1:N
%             A_sur(ii,:) = shrink(B(ii,:),2*lambda/num_may);
%             costo(ii) = norm(A_sur(ii,:),1);
%         end
%         [tau_1] = max(costo);
        [ A_newt, tau_opt,iter] = solve_l1_search_newton( B,lambda, tau_1);
        error(zz,ll) = tau_opt-tau_1;
        ratio(zz,ll) = error(zz,ll)/lambda;    
        mayores(zz,ll) = sum(nb>tau_opt);
    end

end

errores = mean(error);
plot(possible_lambda,errores)