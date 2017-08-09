function [ x ] = xupdate_1inf_overlap( x,z,v,y,rho,H)
%BPDN_XUPDATE Summary of this function goes here
%   Detailed explanation goes here
    lhs = @(w) (w/rho + (H'*H*w));
    b = y/rho + H'*(z-v);
    x = cgs(lhs,b,1E-5,20,[],[],x);

end

