function [ s ] = sigmoid( v, beta,lambda )
%SIGMOID Summary of this function goes here
%   Detailed explanation goes here
        % s(x) = 1 / 1 - exp(-b(x-a))
%     s = 1./(1+exp(-b*(x-a)));
    
    
    s = 0.5*((beta * (v-lambda))./sqrt((1 + (beta^2) * (v-lambda).^2))+1); 
end

