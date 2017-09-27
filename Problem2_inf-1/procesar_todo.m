clc;
clear;
close all;


% N_list = [2000, 5000,10000];
% M_list = [100, 200, 300];

N_list = [5000];
M_list = [200];

for z=1:length(N_list)
N = N_list(z);
M = M_list(z);

disp([ num2str(N) ' x ' num2str(M)])
gen_results_todo(N,M)

end