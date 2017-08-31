clc
clear;
close all;


N = 5000;
M = 1000;
B = (rand(N,M)-0.5);
lambda =0.2;
norma_B = max(sum(abs(B),2));
    [~,ib]=max(sum(abs(B),2));
    tau_0 = norm(shrink(B(ib,:),lambda),1);
tau_range = tau_0:0.01:tau_0*2;

f = @(tau) search_function_l1(tau,B,lambda);
% f_acc = @(tau)search_function_l1(tau,B,lambda)/tau;
% f_acc2 = @(tau) search_function_l1(tau,B,lambda)*tau;
f_expandido = @(tau) f_prueba( tau,B,lambda );
%     
% for ii=1:length(tau_range)
% %     f_normal(ii) = f(tau_range(ii));
% %     acc_f(ii) = f_acc(tau_range(ii));
% %     acc_f2(ii) = f_acc2(tau_range(ii));
%     test_f(ii) = f_expandido(tau_range(ii));
% end
% 
% plot(tau_range,f_normal)
% hold on
% % plot(tau_range,acc_f)
% % plot(tau_range,acc_f2)
% plot(tau_range,test_f)
% % legend('f','f/tau','f*tau')
% 
% legend('f','f prueba')
% 
% sum(abs(test_f-f_normal))

f_newt = @(tau) search_function_l1_der(tau,B,lambda);
f_newt2 = @(tau) search_function_l1_der_2(tau,B,lambda);

max_iter=1000;
tol = 1E-10;
tol_u = 1E-10;


% 
% 
% 

% 
tic
[tau_opt,x_hist,error,iter]=steff_amat(f,tau_0,tol,max_iter,tol_u);
toc

% tic
% [tau_opt2,x_hist2,error2,iter2]=steff_amat(f,0,tol,max_iter,tol_u);
% toc
% error
% error2
% 
% iter
% iter2
tic
[tau_newt,x_hist_newt,error_newt,iter_newt]=root_newt(f_newt,tau_0,tol,max_iter);
toc

% tic
% [tau_newt2,x_hist_newt2,error_newt2,iter_newt2]=root_newt(f_newt,0,tol,max_iter);
% toc

% tic
% [tau_opt_acc,x_hist_acc,error_acc,iter_acc]=steff_amat(f_acc,tau_0,tol,max_iter,tol_u);
% toc
% tic
% [tau_opt_acc2,x_hist_acc2,error_acc2,iter_acc2]=steff_amat(f_acc2,tau_0,tol,max_iter,tol_u);
% toc
% 
error
error_newt
% error_newt2
% 
iter
iter_newt
% iter_newt2
% iter_acc2