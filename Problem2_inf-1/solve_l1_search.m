function [ A_prueba, costo, costo0,tau_opt] = solve_l1_search( B,lambda,possible_t)
%SOLVE_L1_SEARCH Summary of this function goes here
%   Detailed explanation goes here
    norm0 = @(u) sum(abs(u)>0.0001,2);

    jj=1;
    [N,M] = size(B);
    costo = zeros(length(possible_t),1);
    costo0 = zeros(length(possible_t),1);
    for t = possible_t
        [A_prueba] = loop_projL1Mich(B, t/lambda, 20);
        costo(jj) = 0.5*sum_square( B(:) - A_prueba(:)) + lambda*max(sum(abs(A_prueba),2));
        costo0(jj) = 0.5*sum_square( B(:) - A_prueba(:)) + lambda*max(norm0(A_prueba));
        a_store{jj} = A_prueba;
        jj = jj+1;
    end

[tau_opt,ind] = min(costo);
A_prueba = a_store{ind};

end

