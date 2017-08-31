clc;
clear;
close all;

N = 500;
M = 200;
% rng(3);

shrink = @(u,ll) sign(u).*max(abs(u)-ll,0);

num_real = 100;

errores = zeros(num_real,2);
tiempo = zeros(num_real,2);

for zz=1:num_real
%     rng(5*zz)
    zz
    
    B = (rand(N,M)-0.5);
    lambda =0.1;

    tic
    norma_B = max(sum(abs(B),2));
    [~,ib]=max(sum(abs(B),2));
    tau_1 = norm(shrink(B(ib,:),lambda),1);
    [ A_stef, tau_opt_stef,iter] = solve_l1_search_steffenson( B,lambda, tau_1);
    X_stef = B-A_stef;
    tiempo(zz,2) = toc;
    errores(zz,2) = abs(compute_mixed_norm(X_stef,1,inf)-lambda);
    iter_num(zz,2) = iter;
    clear X_stef A_stef;
    
    tic
    norma_B = max(sum(abs(B),2));
    [~,ib]=max(sum(abs(B),2));
    tau_1 = norm(shrink(B(ib,:),lambda),1);
    [ A_newt, tau_opt_stef,iter] = solve_l1_search_newton( B,lambda, tau_1);
    X_newt = B-A_newt;
    tiempo(zz,3) = toc;
    errores(zz,3) = abs(compute_mixed_norm(X_newt,1,inf)-lambda);
    iter_num(zz,3) = iter;
    clear X_newt A_newt;
    
    tic
    [ X_sra, theta_opt, iter ] = solve_sra(B,lambda);
    tiempo(zz,1) = toc;
    errores(zz,1) =abs(compute_mixed_norm(X_sra,1,inf)-lambda);
    iter_num(zz,1) = iter;
    clear X_sra

end