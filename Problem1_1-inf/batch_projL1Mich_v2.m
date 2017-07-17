function[x, l, loops] = batch_projL1Mich_v2(B, tau, nMaxIter, myErr)

if nargin < 4
  myErr = 1e-8;
end
[nrows,ncols] = size(B);


B = B';

b = B(:);
% L = kron(speye(nrows),ones(1,ncols));
% L = sparse(kron(eye(nrows),ones(1,ncols)));
sum_blocks1 = @(M) sum(reshape(M,[length(M)/nrows,nrows]))';

% sum_blocks = @(M) L*M;

s0 = sign(b);
bnorm = sum_blocks1(s0.*b); % vector of s(j)'*b(j) values
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
      sN = sum_blocks1(s.*s);
      sb = sum_blocks1(s.*b);
      
% =================================

for k=1:nMaxIter
    
    indic = (sN == sN0);
    
    if all(indic)
      break;
    end    
    
    l = ( (sb - tau)./(sN) );
          
    sN0 = sN;
    
    s = s0.*(bAbs > kron(l,ones(ncols,1)));
    sN(~indic) = sum_blocks(s.*s,~indic,ncols,nrows);
    sb(~indic) = sum_blocks(s.*b,~indic,ncols,nrows);
     
end

x = s.*max(0, (bAbs - kron(l,ones(ncols,1))));  
x = reshape(x,[ncols,nrows])';
loops = k;
end

function s = sum_blocks(M,indices,ncols,nrows)
    aa = reshape(M,[ncols,nrows]);
    s = sum(aa(:,indices))';
end



    