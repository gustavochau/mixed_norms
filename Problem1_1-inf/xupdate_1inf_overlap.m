function [ x ] = xupdate_1inf_overlap( x,z,A,B,c,v,rho,params, AT )
%BPDN_XUPDATE Summary of this function goes here
%   Detailed explanation goes here
    y = params.y;
    mat_form = params.mat_form;
    G_aux = params.G_aux;
    
    z = vec(mat_form(z)');
    v = vec(mat_form(v)');
    
    lhs = @(w) (w + rho*(G_aux'*G_aux*w));
    b = y + rho*G_aux'*(z+v);
    x = cgs(lhs,b,1E-5,20,[],[],x);

end

