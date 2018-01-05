function [ gprime ] = der_function(lambda, v, beta, tau)
%DER_FUNCTION Summary of this function goes here
%   Detailed explanation goes here

%     gprime = sum((lambda-v).*sigmoid(v,beta,lambda)) + 0.5*beta*sum(((v-lambda).^2).* sigmoid(v,beta,lambda).* (1-sigmoid(v,beta,lambda)))+tau;
    %gprime = gprime -tau + sum((v-2*lambda).*sigmoid(v,beta,lambda)) + beta*sum((v*lambda-lambda^2).*sigmoid(v,beta,lambda).* (1-sigmoid(v,beta,lambda)));
    der_sig = -(0.5*beta*sqrt(beta^2*(v - lambda).^2 + 1))./((beta^2 * (v-lambda).^2+1 ).^2);
    gprime = sum(-(v-lambda).*sigmoid(v,beta,lambda)) + 0.5*beta*sum(((v-lambda).^2).* der_sig)+tau;
end

