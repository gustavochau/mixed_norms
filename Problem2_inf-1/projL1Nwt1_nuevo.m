function[x, l, loops] = projL1Nwt1_nuevo(b, tau, nMaxIter, myErr)
loops=1;

if nargin < 4
  myErr = 1e-8;
end
if tau==0
    x = b*0;
    l = 0;
    return;
end
l=0;
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
    
    alpha = ( (sb - tau)/(sN) )/l0;
          
    l = l0*(1.0/(2-alpha));
        
    sN0 = sN;
    
    s = s0.*(bAbs>l);
    sN = s'*s;
        
    if sN == 0
      l = ( (sb - tau)/(sN0) );
      s = s0.*(bAbs>l);
      sN = s'*s;      
    end
    
    sb = s'*b;                             
    l0 = l;
      
end

l = ( (sb - tau)/(sN) );
    
x = s.*max(0, bAbs - l);  
loops = k;

    