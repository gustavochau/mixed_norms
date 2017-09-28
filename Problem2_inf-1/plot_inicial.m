clc;
clear;
close all;
load('p_inciial.mat')
gamma =[0.0001:0.0001:0.001];
taus = squeeze(mean(taus,1))';

plot(gamma,taus,'linewidth',1.5)
legend('tau_0','tau^*')
xlabel('gamma')
ylabel('Value')

iter = squeeze(resumen(:,2,:))';
figure
plot(gamma,iter,'linewidth',1.5)
xlabel('gamma')
ylabel('number of iterations')
legend('Starting from 0','Starting from tau_0')
