clc;
clear;
close all;

tau = 0.5;
b = randn(50,1);
[x, l] =    projL1Mich(b, tau,20);
v=abs(b);
beta=100;

domain = 0:0.001:4*l;
for ii=1:length(domain)
    lambda=domain(ii);
    gprime(ii) = der_function(lambda, v, beta, tau);
end
fun = @(w) der_function(w, v, beta, tau);
l2 = fzero(fun,0)

plot(domain,gprime)
hold on
scatter(l,0,'Marker','x')


syms lambda v tau beta


eqn = -(v-lambda)*0.5*((beta * (v-lambda))/sqrt((1 + (beta^2) * (v-lambda)^2))+1) + 0.5*beta*((v-lambda)^2)* -(0.5*beta*sqrt(beta^2*(v - lambda)^2 + 1))/((beta^2 * (v-lambda)^2+1 )^2)+tau == 0;
% eqn = (lambda-v)*(1/(1+exp(beta*(v-lambda)))) +0.5*beta*((v-lambda)^2)*(1/(1+exp(beta*(v-lambda))))*(1-(1/(1+exp(beta*(v-lambda)))))+tau == 0;

%  0.5*(  (\beta * x)/(1 + \beta^2 * x^2) - 1)

%    gprime = sum((lambda-v).*sigmoid(v,beta,lambda)) + 0.5*beta*sum(((v-lambda).^2).* sigmoid(v,beta,lambda).* (1-sigmoid(v,beta,lambda)))+tau;


 solx = solve(eqn, lambda);