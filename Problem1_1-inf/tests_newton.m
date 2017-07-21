clc;
clear;
close all;

% 
% % square
% size_matrices{1} = [50 50]; 
% size_matrices{2} = [100 100];
% size_matrices{3} = [500 500];
% 
% % tall
% size_matrices{4} = [200 100];
% size_matrices{5} = [500 100];
% size_matrices{6} = [1000 100];
size_matrices{1} = [5000 100];

% wide
% size_matrices{8} = [100 200];
% size_matrices{9} = [100 500];
% size_matrices{10} = [100 1000];
size_matrices{2} = [100 5000];

metrics = zeros(length(size_matrices),3);

for zz = 1:length(metrics)
    disp(['case' num2str(zz)])
    sm = size_matrices{zz};
    
    N = sm(1);
    M = sm(2);
% N = 100;
% N2 = 100;
B = rand(N,M);
lambda = 0.1;
num_rep = 1000;
loop_time=zeros(num_rep,1);
block_time=zeros(num_rep,1);
differ=zeros(num_rep,1);
for ii=1:num_rep
tic
[ A, l1, loops ] = loop_projL1AccNewton( B, lambda, 20 );
loop_time(ii) = toc;

tic
[ A2, l2, loops2 ] = batch_projL1AccNewton( B, lambda, 20 );
block_time(ii) = toc;

% 
% tic
% [ A3, l3, loops3 ] = batch_projL1Mich_v2( B, lambda, 20 );
% block_time2(ii) = toc;
differ(ii) = norm(A2-A,'fro');

end
metrics(zz,1) = mean(loop_time);
metrics(zz,2) = mean(block_time);
metrics(zz,3) = max(abs(differ));
% max(differ)
% disp(['loops: ' num2str(mean(loop_time))])
% disp(['block: ' num2str(mean(block_time))])
% disp(['block v2: ' num2str(mean(block_time2))])


end
