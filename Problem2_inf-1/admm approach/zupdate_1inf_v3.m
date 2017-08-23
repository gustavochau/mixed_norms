function [ z ] = zupdate_1inf_v3( alpha,v,rho,lambda,LS)
%BPDN_XUPDATE Summary of this function goes here
%   Detailed explanation goes here
    b = LS(alpha) + v;
%     z = b-projL1AccNewton(b, lambda/rho, 40);
    z = b - projL1Mich(b, lambda/rho, 40);
end

