clc;
clear;
close all;

N = 10000;
realizations=1000;
rng(3)



for zz=1:realizations
    rng(zz) % for repeatibility
    b = rand(N,1)-0.5; %uniform between -0.5 and 0.5
    tau = 0.1;
    
%     % cvx solution
%     cvx_begin 
%         variable a_cvx(N,1)
%         minimize( 0.5*sum_square( b - a_cvx))
%         subject to
%             norm(a_cvx,1) <= tau
%     cvx_end


    
    tic
    [a_mich,lam_mich,l_mich] = projL1Mich(b, tau, 20); %michelot
    tiempo(zz,1) = toc;
    tic
    [a_newt,lam_newt,l_newt] = projL1Nwt1(b, tau, 20); % newton r fixed
    tiempo(zz,2)=toc;
    tic
    [a_newt1,lam_newt1,l_newt1] = projL1AccNewton_final(b, tau, 20); % newton v2
    tiempo(zz,3)=toc;
    
%     % Error wrt cvx
%     error(zz,1) = max(abs(a_mich - a_cvx));
%     error(zz,2) = max(abs(a_newt - a_cvx));
%     error(zz,3) = max(abs(a_newt1 - a_cvx));

    % constraint violation ||a||_1 - tau
%     const_viol(zz,1) = abs(norm(a_cvx,1)-tau);
    const_viol(zz,2) = abs(norm(a_mich,1)-tau);
    const_viol(zz,3) = abs(norm(a_newt,1)-tau);
    const_viol(zz,4) = abs(norm(a_newt1,1)-tau);
    
    
end