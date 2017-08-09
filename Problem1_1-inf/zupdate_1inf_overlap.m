function [ z ] = zupdate_1inf_overlap( x,v,rho,lambda,G,mat_form)
%UPDATE_C_L1BALL Summary of this function goes here
%   Detailed explanation goes here
    X = diag(x);
    V = mat_form(v);
    
    target =  G*X+V;
    Z =  solve_l1_inf_prox_loop(target, lambda/rho, 20);
    
    z = vec(Z');

end

