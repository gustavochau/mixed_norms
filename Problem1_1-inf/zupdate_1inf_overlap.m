function [ z ] = zupdate_1inf_overlap( x,z,A,B,c,v,rho,params, AT )
%UPDATE_C_L1BALL Summary of this function goes here
%   Detailed explanation goes here
    mat_form = params.mat_form;    
    target =  A(diag(x))-mat_form(v);
    Z =  solve_l1_inf_prox_loop(target, params.lambda/rho, 20);
    z = Z(:);

end

