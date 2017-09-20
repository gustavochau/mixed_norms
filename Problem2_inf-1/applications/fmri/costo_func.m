function [ costo ] = costo_func( U,X,Y )
%COSTO_FUNC Summary of this function goes here
%   Detailed explanation goes here
    costo = 0;
    num_tasks = size(Y,2);
    
    for ii=1:num_tasks
        costo = 0.5*norm(Y(:,ii)-X*U(:,ii));
    end

end

