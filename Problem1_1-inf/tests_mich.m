clc;
clear;
close all;

N = 10000;
N2 = 200;
B = rand(N,N2);
lambda = 0.1;

for ii=1:1000
tic
[ A, l1, loops ] = loop_projL1Mich( B, lambda, 20 );
loop_time(ii) = toc;
% error_cvx = norm(A-A_cvx,'fro')/N/N

tic
[ A2, l2, loops2 ] = batch_projL1Mich( B, lambda, 20 );
block_time(ii) = toc;

differ(ii) = norm(A-A2,'fro');

end

max(differ)
disp(['loops: ' num2str(mean(loop_time))])
disp(['block: ' num2str(mean(block_time))])

