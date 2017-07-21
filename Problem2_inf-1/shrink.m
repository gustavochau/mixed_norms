function [ shrinked_x ] = shrink( x, thresh )
%SHRINK Aplica el operador shrinkage usando el umbral thresh
    shrinked_x = sign(x).*max(abs(x)-thresh,0);
end

