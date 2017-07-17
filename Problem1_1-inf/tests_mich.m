clc;
clear;
close all;

N = 500;
N2 = 100;
B = rand(N,N2);
lambda = 0.1;

for ii=1:1
% tic
% [ A, l1, loops ] = loop_projL1Mich( B, lambda, 20 );
% loop_time(ii) = toc;
% % error_cvx = norm(A-A_cvx,'fro')/N/N

tic
[ A2, l2, loops2 ] = batch_projL1Mich( B, lambda, 20 );
block_time(ii) = toc;

% 
tic
[ A3, l3, loops3 ] = batch_projL1Mich_v2( B, lambda, 20 );
block_time2(ii) = toc;
differ(ii) = norm(A2-A3,'fro');

end

max(differ)
% disp(['loops: ' num2str(mean(loop_time))])
disp(['block: ' num2str(mean(block_time))])
disp(['block v2: ' num2str(mean(block_time2))])

