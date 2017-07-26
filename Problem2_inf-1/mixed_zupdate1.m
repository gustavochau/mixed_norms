function [ z ] = mixed_zupdate1( x,z,A,B,c,v,rho,params, AT )
%UPDATE_C_L1BALL Summary of this function goes here
%   Detailed explanation goes here
    L = params.L;
    target = L*abs(x);
    z =  target- batch_projL1AccNewton(target, 2*params.lambda/rho, 20);


end

