function [ A, l, loops ] = loop_projL1AccNewt( B, lambda, nMaxIter )

[nrows,ncols] = size(B);
l = zeros(nrows,1);

loops = zeros(nrows,1);
A = zeros(nrows,ncols);

for ii=1:nrows
    [x, loops(ii),l(ii)] = projL1Nwt1_nuevo(B(ii,:)', lambda, nMaxIter);
    A(ii,:) = x';
end

end

