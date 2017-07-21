function [ z ] = update_c_l1ball( alpha,z,A,B,c,v,rho,paramsZ )
%UPDATE_C_L1BALL Summary of this function goes here
%   Detailed explanation goes here
    z = batch_projL1AccNewton(A(alpha), 2*paramsZ.lambda/rho, 20, myErr);


end

