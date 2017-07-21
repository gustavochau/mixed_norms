function [ A, l, loops ] = loop_projL1AccNewton( B, lambda, nMaxIter )

[nrows,ncols] = size(B);
l = zeros(nrows,1);

loops = zeros(nrows,1);
A = zeros(nrows,ncols);

for ii=1:nrows
    [x, l(ii), loops(ii)] = projL1AccNewton(B(ii,:)', lambda, nMaxIter);
    A(ii,:) = x';
end

end

