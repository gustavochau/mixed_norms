function [ x ] =  xupdate_1inf_v3( alpha,z,v,beta,rho,A_list,N,M)
%BPDN_XUPDATE Summary of this function goes here
%   Detailed explanation goes here
%     lhs = @(w) (w + rho*(LS_t*(LS*w)));

    x = zeros(size(beta));
    for ii=1:N
        b = beta( ((ii-1)*M+1): (ii)*M ) + rho* sign(beta( ((ii-1)*M+1): (ii)*M ))*(z(ii)-v(ii));
        x(((ii-1)*M+1):ii*M) = A_list(:,:,ii)'\(A_list(:,:,ii)\b);
    end
    
%     evalc('x = cgs(lhs,b,1E-5,20,[],[],x);');
    
end

