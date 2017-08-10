function [ z ] = zupdate_1inf( alpha,v,rho,lambda,L,S)
%BPDN_XUPDATE Summary of this function goes here
%   Detailed explanation goes here
    b = L*S*alpha + v;
%     z = b-projL1AccNewton(b, lambda/rho, 40);
    z = b - projL1Mich(b, lambda/rho, 40);
end

