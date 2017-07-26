function [ x ] = mixed_xupdate1( x,z,A,B,c,v,rho,params, AT )
%BPDN_XUPDATE Summary of this function goes here
%   Detailed explanation goes here
    beta = params.beta;
    L = params.L;
    for ii=1:5
        A = @(w) L*(sign(x).*w);
        AT = @(w) (L*diag(sign(x)))'*w;
        lhs = @(w) rho*AT(A(w)) + w; % handle al lhs de la ecuación
        b = beta+ rho * AT(z-v ); % Lado derecho de la ecuación
        x = cgs(lhs,b,1E-5,1,[],[],x);
    end
end

