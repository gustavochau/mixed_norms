clc;
clear;
close all;
load('results_pinicial.mat')
gamma =[0.0001:0.0001:0.001];
taus = squeeze(mean(taus,1))';

figure
ha=tight_subplot(1,2,0.1,[ 0.15 0.08],0.05);
axes(ha(1))
plot(gamma,taus,'linewidth',1.5)
legend('\gamma_0','\gamma^*')
xlabel('\alpha')
ylabel('Value')
axis tight
title('(a)')

axes(ha(2))
iter = squeeze(resumen(:,2,:))';
plot(gamma,iter,'linewidth',1.5)
xlabel('\alpha')
ylabel('number of iterations')
legend('Starting from 0','Starting from \gamma_0')
axis tight
title('(b)')