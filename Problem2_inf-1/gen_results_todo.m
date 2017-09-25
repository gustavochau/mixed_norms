function gen_results_todo (N,M)
% 
% N = 5000;
% M = 200;
% rng(3);



num_real = 100;

errores = zeros(num_real,2);
tiempo = zeros(num_real,2);
pp=1;

nzr = @(U) sum(max(abs(U),[],2)> 1E-4);

for  gamma =[0.0001,0.0005,0.001]; %lambda =[0.05,0.1,0.2];
%     gamma
    for zz=1:num_real
        %     rng(5*zz)
        disp(num2str(zz))
        rng(7*zz)
        B = (rand(N,M)-0.5);
         lambda = gamma*compute_mixed_norm(B,1,inf);
%         lambda = gamma*10;
        
        % sra
        tic
        [ X_sra, theta_opt, iter ] = solve_sra(B,lambda);
        tiempo(zz,1) = toc;
        errores(zz,1) =abs(compute_mixed_norm(X_sra,1,inf)-lambda);
        iter_num(zz,1) = iter;
        nonzero(zz,1) = nzr(X_sra)*100/N;
        clear X_sra
        
        % steffensen
        tic
        [ X_stef, theta_opt, iter ] = proj_steffensen(B,lambda);
        tiempo(zz,2) = toc;
        errores(zz,2) =abs(compute_mixed_norm(X_stef,1,inf)-lambda);
        iter_num(zz,2) = iter;
        nonzero(zz,2) = nzr(X_stef)*100/N;
        clear X_stef
        
        % steffensen + pruning
        tic
        [ X_stef_prun, theta_opt, iter ] = proj_steffensen_pruned(B,lambda);
        tiempo(zz,3) = toc;
        errores(zz,3) =abs(compute_mixed_norm(X_stef_prun,1,inf)-lambda);
        iter_num(zz,3) = iter;
        nonzero(zz,3) = nzr(X_stef_prun)*100/N;
        clear X_stef_prun
        
        % newton
        tic
        [X_newt, theta_opt, iter ] = proj_newton( B, lambda);
        tiempo(zz,4) = toc;
        errores(zz,4) = abs(compute_mixed_norm(X_newt,1,inf)-lambda);
        iter_num(zz,4) = iter;
        nonzero(zz,4) = nzr(X_newt)*100/N;
        clear X_newt;
        
        % newton + pruning
        tic
        [X_newt_prun, theta_opt, iter ] = proj_newton_pruned( B, lambda);
        tiempo(zz,5) = toc;
        errores(zz,5) = abs(compute_mixed_norm(X_newt_prun,1,inf)-lambda);
        iter_num(zz,5) = iter;
        nonzero(zz,5) = nzr(X_newt_prun)*100/N;
        clear X_newt_prun;
        
    end
    err_cell{pp} = errores;
    iter_cell{pp} = iter_num;
    tiempo_cell{pp} = tiempo;
    nz_cell{pp} = nonzero;
    resumen(:,:,pp) = [mean(errores)' mean(iter_num)' mean(tiempo)' mean(nonzero)'];
    pp=pp+1;
end
save(['results_' num2str(N) 'x' num2str(M) '_lambda' num2str(lambda) '.mat'],'resumen','err_cell','iter_cell','tiempo_cell','nz_cell')

end