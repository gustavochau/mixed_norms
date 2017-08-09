function [ mix_norm ] = compute_mixed_norm( X, norm_1, norm_2 )
%COMPUTE_MIXED_NORM Summary of this function goes here
%   Detailed explanation goes here
    [nrows,ncols] = size(X);
    norm_rows = zeros(nrows,1);
    for ii=1:nrows
        norm_rows(ii) = norm(X(ii,:),norm_2);
    end
    mix_norm = norm(norm_rows,norm_1);

end

