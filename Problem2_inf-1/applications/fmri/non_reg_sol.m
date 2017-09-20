function [ W ] = non_reg_sol( X,Y )
%NON_REG_SOL Summary of this function goes here
%   Detailed explanation goes here
    num_tasks = size(Y,2);
    num_feats = size(X,2);
    W = zeros(num_feats,num_tasks);
    opA = @(U)X'*(X*U);
    for ii=1:num_tasks
         W(:,ii) = cgs(opA,X'*Y(:,ii),1E-5,2000);
    end


end

