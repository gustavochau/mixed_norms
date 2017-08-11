function[x, l, loops] = projL1AccNewton(b, tau, nMaxIter, myErr)

if nargin < 4
  myErr = 1e-8;
end

s0 = sign(b);
bnorm = s0'*b;
bAbs = s0.*b;


sN0 = length(b(:));


xnorm = bnorm;


if bnorm <= tau % if b belongs to ball then b is the solution
  x = b;
  return;
end

% ==================================================
%  Init

      l0 = (bnorm - tau)/sN0;      
            
      s = s0.*(bAbs>l0);
      sN = s'*s;
      sb = s'*b;
      
      
% =================================

r = 1;

for k=1:nMaxIter

      
    if sN == sN0
      break;
    end    
    
    alpha = ( (sb - tau)/(sN) )/l0;

    if( k >= 2)
      r = 0.5/(alpha - 1.0);
      r = min( [r, 6] );
      r = max( [r, 1] );
    end
    
    l = l0*(1 - ((r-1)/r)*alpha) / ( (r+1)/r - alpha);
    
    
    sN0 = sN;
    
    s = s0.*(bAbs>l);
    sN = s'*s;
    sb = s'*b;
                             
    l0 = l;
      
end

    
x = s.*max(0, bAbs - l);  
loops = k;

    