function[x, l, loops] = batch_projL1Mich(B, tau, nMaxIter, myErr)

if nargin < 4
  myErr = 1e-8;
end
[nrows,ncols] = size(B);


B = B';

b = B(:);
sum_blocks = @(M) sum(reshape(M,[ncols,nrows]))';

s0 = sign(b);
bnorm = sum_blocks(s0.*b); % vector of s(j)'*b(j) values
bAbs = s0.*b;


sN0 = ones(nrows,1)*ncols;

xnorm = bnorm;


if bnorm <= tau
  x = reshape(b,[ncols,nrows])';
  return;
end

% ==================================================
%  Init

      l0 = (bnorm - tau)./sN0;      
            
      s = s0.*( bAbs > kron(l0,ones(ncols,1)));
      sN = sum_blocks(s.*s);
      sb = sum_blocks(s.*b);
      
% =================================

for k=1:nMaxIter
    
    if all(sN == sN0)
      break;
    end    
    
    l = ( (sb - tau)./(sN) );
          
    sN0 = sN;
    
    s = s0.*(bAbs > kron(l,ones(ncols,1)));
    sN = sum_blocks(s.*s);
    sb = sum_blocks(s.*b);
     
end

x = s.*max(0, (bAbs - kron(l,ones(ncols,1))));  
x = reshape(x,[ncols,nrows])';
loops = k;

    