function [ x ] = xupdate_1inf_v2( x,z,v,beta,rho,A,A_t,LS_t)
%BPDN_XUPDATE Summary of this function goes here
%   Detailed explanation goes here
%     lhs = @(w) (w + rho*(LS_t*(LS*w)));
    b = beta + rho*LS_t*(z-v);
%     evalc('x = cgs(lhs,b,1E-5,20,[],[],x);');
    x = A_t\(A\b);
end

