function [ x ] = xupdate_1inf( x,z,v,beta,rho,L,S)
%BPDN_XUPDATE Summary of this function goes here
%   Detailed explanation goes here
    lhs = @(w) (w + rho*(S'*(L'*(L*(S*w)))));
    b = beta + rho*S'*(L'*(z-v));
    evalc('x = cgs(lhs,b,1E-5,20,[],[],x);');

end

