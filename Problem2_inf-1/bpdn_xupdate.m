function [ x ] = bpdn_xupdate( x,z,A,B,c,v,rho,params )
%BPDN_XUPDATE Summary of this function goes here
%   Detailed explanation goes here
    Phi = params.K;
    Phi_T = params.KT;
    In = params.In;
   
    lhs = @(w) Phi_T(Phi(w)) + rho*w; % handle al lhs de la ecuación
    b = Phi_T(In)-rho*(v-z); % Lado derecho de la ecuación
    x = cgs(lhs,b,1E-5,30,[],[],x);
    
end

