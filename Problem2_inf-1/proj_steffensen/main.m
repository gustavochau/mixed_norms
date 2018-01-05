clc;
clear;
close all;


% size of matrix
N = 50;
M = 20;

num_real = 100; % number of realizations

errores = zeros(num_real,2); % constraint violation
tiempo = zeros(num_real,2); % run time
iter_num = zeros(num_real,2); % number of iteraitons
nonzero= zeros(num_real,2); %number of nonzero rows
pp=1;

nzr = @(U) sum(max(abs(U),[],2)> 1E-4); % anonymoous function for computing number of nonzero rows

for  gamma =[0.0001,0.0005,0.001]; %lambda =[0.05,0.1,0.2];
    
    for zz=1:num_real
        disp(num2str(zz))
        rng(7*zz)
        B = (rand(N,M)-0.5); %generate matrix
        tau = gamma*compute_mixed_norm(B,1,inf); % tau depending on mixed norm of B
       
        % sra
        tic
        [ X_sra, theta_opt, iter ] = solve_sra(B,tau);
        tiempo(zz,1) = toc;
        errores(zz,1) =abs(compute_mixed_norm(X_sra,1,inf)-tau);
        iter_num(zz,1) = iter;
        nonzero(zz,1) = nzr(X_sra)*100/N;
        clear X_sra
     
        % steffensen + pruning
        tic
        [ X_stef_prun, theta_opt, iter ] = proj_steffensen_pruned(B,tau);
        tiempo(zz,2) = toc;
        errores(zz,2) =abs(compute_mixed_norm(X_stef_prun,1,inf)-tau);
        iter_num(zz,2) = iter;
        nonzero(zz,2) = nzr(X_stef_prun)*100/N;
        clear X_stef_prun
                
    end
    err_cell{pp} = errores;
    iter_cell{pp} = iter_num;
    tiempo_cell{pp} = tiempo;
    nz_cell{pp} = nonzero;
    resumen(:,:,pp) = [mean(errores)' mean(iter_num)' mean(tiempo)' mean(nonzero)']; % average over all realizations
    pp=pp+1;
end