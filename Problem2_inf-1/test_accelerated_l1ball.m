clc
clear;
close all;


N = 1000;
u = (rand(N,1)-0.5);
tau =0.1;

l_range = 0.01:0.01:2;

f = @(l) sum(max(abs(u)-l,0))-tau;
f2 = @(l) l*f(l);
    
for ii=1:length(l_range)
    f_normal(ii) = f(l_range(ii));
    acc_f(ii) = f2(l_range(ii));
end

plot(l_range,f_normal)
hold on
plot(l_range,acc_f)
% plot(l_range,acc_f2)

legend('f','f*tau','f/tau')

tau_0=0.37;
max_iter=1000;
tol_u = 1E-10;
tol = 1E-10;
tic
[tau_opt,x_hist,error,iter]=steff_amat(f,tau_0,tol,max_iter,tol_u);
toc
tic
[tau_opt_acc,x_hist_acc,error_acc,iter_acc]=steff_amat(f2,tau_0,tol,max_iter,tol_u);
toc
% tic
% [tau_opt_acc2,x_hist_acc2,error_acc2,iter_acc2]=steff_amat(f_acc2,tau_0,tol,max_iter,tol_u);
% toc

iter
iter_acc
% iter_acc2