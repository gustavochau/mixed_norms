clc;
clear;
close all;

N = 100;
M = 100;



num_test = 1;

% for pp=1:num_test
% rng(3);
B = (rand(N,M));
% B = ones(N,M);

% B= [1  1 1; 0.8 0.8 0.8; 0.8 0.8 0.8]
% B(1,:) = B(1,:)*2;
kk=1;
for lambda = 0.01:0.01:3


shrink = @(u,ll) sign(u).*max(abs(u)-ll,0);

for ii=1:N
    A_sur(ii,:) = shrink(B(ii,:),lambda);
    costo(ii) = norm(A_sur(ii,:),1);
end
[tau,ind] = max(costo);


[~,jj] = max(sum(abs(B),2));

tau2 = norm(shrink(B(jj,:),lambda),1);
error(kk) = tau-tau2;
disp(error)
kk=kk+1;
end
plot(0.01:0.01:3,error)