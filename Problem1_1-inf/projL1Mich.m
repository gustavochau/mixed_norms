function[x, l,loops] = projL1Mich(b, tau, nMaxIter, myErr)
loops=1;
if nargin < 4
  myErr = 1e-8;
end

s0 = sign(b);
bnorm = s0'*b;
bAbs = s0.*b;


sN0 = length(b(:));


xnorm = bnorm;


if bnorm <= tau
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

for k=1:nMaxIter

      
    if sN == sN0
      break;
    end    
    
    l = ( (sb - tau)/(sN) );
          
    sN0 = sN;
    
    s = s0.*(bAbs>l);
    sN = s'*s;
    sb = s'*b;
                             
      
end

    
x = s.*max(0, bAbs - l);  
loops = k;

    