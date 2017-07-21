function [ z ] = bpdn_zupdate( x,z,A,B,c,v,rho,params )
%BPDN_XUPDATE Summary of this function goes here
%   Detailed explanation goes here
      lambda = params.lambda;  
      z =  shrink(x+v, lambda/rho);
end

