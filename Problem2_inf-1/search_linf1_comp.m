clc;
clear;
close all;




N = 5000;
M = 200;
% rng(3);

shrink = @(u,ll) sign(u).*max(abs(u)-ll,0);

num_real = 5;

errores = zeros(num_real,2);
tiempo = zeros(num_real,2);
pp=1;
for lambda =[0.05,0.1,0.2];
    lambda
    for zz=1:num_real
        %     rng(5*zz)
        disp(num2str(zz))
        
        B = (rand(N,M)-0.5);
        
        %     tic
        %     norma_B = max(sum(abs(B),2));
        %     [~,ib]=max(sum(abs(B),2));
        %     tau_1 = norm(shrink(B(ib,:),lambda),1);
        %     [ A_stef, tau_opt_stef,iter] = solve_l1_search_steffenson( B,lambda, tau_1);
        %     X_stef = B-A_stef;
        %     tiempo(zz,2) = toc;
        %     errores(zz,2) = abs(compute_mixed_norm(X_stef,1,inf)-lambda);
        %     iter_num(zz,2) = iter;
        %     clear X_stef A_stef;
                
        tic
        [ X_sra, theta_opt, iter ] = solve_sra(B,lambda);
        tiempo(zz,1) = toc;
        errores(zz,1) =abs(compute_mixed_norm(X_sra,1,inf)-lambda);
        iter_num(zz,1) = iter;
        clear X_sra
        
        tic
        for ii=1:N
            A_sur(ii,:) = shrink(B(ii,:),lambda);
            costo(ii) = norm(A_sur(ii,:),1);
        end
        [tau_1] = max(costo);
        [ A_newt, tau_opt_stef,iter] = solve_l1_search_newton( B,lambda, tau_1);
        X_newt = B-A_newt;
        tiempo(zz,2) = toc;
        errores(zz,2) = abs(compute_mixed_norm(X_newt,1,inf)-lambda);
        iter_num(zz,2) = iter;
        clear X_newt A_newt;
        
        tic
        for ii=1:N
            A_sur(ii,:) = shrink(B(ii,:),lambda);
            costo(ii) = norm(A_sur(ii,:),1);
        end
        [tau_1] = max(costo);
        [ A_newt_prun, tau_opt_stef,iter] = solve_l1_search_newton_pruned( B,lambda, tau_1);
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