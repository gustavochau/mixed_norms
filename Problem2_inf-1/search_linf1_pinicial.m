clc;
clear;
close all;




N = 2000;
M = 100;
% rng(3);

num_real = 100;

errores = zeros(num_real,2);
tiempo = zeros(num_real,2);
pp=1;

nzr = @(U) sum(max(abs(U),[],2)> 1E-4);

for  gamma =[0.0001:0.0001:0.001]; %lambda =[0.05,0.1,0.2];
%     gamma
    for zz=1:num_real
        %     rng(5*zz)
        disp(num2str(zz))
        
        B = (rand(N,M)-0.5);
        lambda = gamma*compute_mixed_norm(B,1,inf);
        
        % steffensen inital 0
        tic
        [ X_stef, tau_opt, iter ] = proj_steffensen(B,lambda,0);
        tiempo(zz,1) = toc;
        errores(zz,1) =abs(compute_mixed_norm(X_stef,1,inf)-lambda);
        iter_num(zz,1) = iter;
        nonzero(zz,1) = nzr(X_stef)*100/N;
        clear X_stef
 
        % steffensen initial point
        tic
        [ X_stef, tau_opt, iter, tau_1 ] = proj_steffensen(B,lambda);
        tiempo(zz,2) = toc;
        errores(zz,2) =abs(compute_mixed_norm(X_stef,1,inf)-lambda);
        iter_num(zz,2) = iter;
        nonzero(zz,2) = nzr(X_stef)*100/N;
        tau_1_hist(zz,1) = tau_1;
        tau_opt_hist(zz,1) = tau_opt;
        clear X_stef
    end
%     err_cell{pp} = errores;
%     iter_cell{pp} = iter_num;
%     tiempo_cell{pp} = tiempo;
    taus(:,:,pp) = [tau_1_hist, tau_opt_hist];
    resumen(:,:,pp) = [mean(errores)' mean(iter_num)' mean(tiempo)' mean(nonzero)'];
    pp=pp+1;
end
% save(['results_' num2str(N) 'x' num2str(M) '_lambda' num2str(lambda) '.mat'],'resumen','err_cell','iter_cell','tiempo_cell','nz_cell')
