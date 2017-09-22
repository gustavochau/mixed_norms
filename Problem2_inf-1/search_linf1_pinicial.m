clc;
clear;
close all;




N = 500;
M = 100;
% rng(3);

shrink = @(u,ll) sign(u).*max(abs(u)-ll,0);

num_real = 100;

errores = zeros(num_real,2);
tiempo = zeros(num_real,2);
pp=1;
for lambda =[0.05,0.1,0.2];
    lambda
    for zz=1:num_real
        %     rng(5*zz)
        disp(num2str(zz))
        
        B = (rand(N,M)-0.5);
        
              
        tic
        for ii=1:N
            A_sur(ii,:) = shrink(B(ii,:),lambda);
            costo(ii) = norm(A_sur(ii,:),1);
        end
        [tau_1] = max(costo);
        [ A_newt_prun, tau_opt_stef,iter] = solve_l1_search_newton_pruned( B,lambda, 0);
        X_newt_prun = B-A_newt_prun;
        tiempo(zz,3) = toc;
        errores(zz,3) = abs(compute_mixed_norm(X_newt_prun,1,inf)-lambda);
        iter_num(zz,3) = iter;
        clear X_newt_prun A_newt_prun;
        
    end
    resumen(:,:,pp) = [mean(errores)' mean(iter_num)' mean(tiempo)'];
    pp=pp+1;
    save(['results_' num2str(N) 'x' num2str(M) '_lambda' num2str(lambda) '.mat'],'iter_num','errores','tiempo','resumen')

end