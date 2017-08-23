function [ A_prueba, costo ] = solve_l0_search( B,lambda,possible_t)
%SOLVE_L0_SEARCH Summary of this function goes here
%   Detailed explanation goes here
    jj=1;
    [N,M] = size(B);
    A_prueba = zeros(N,M);
    costo = zeros(N,1);
    for t = possible_t
        for ii=1:N
            [A_prueba(ii,:)] = MP_l0(B(ii,:), t/lambda);
        end
        costo(jj) = 0.5*sum_square( B(:) - A_prueba(:)) + lambda*max(sum(abs(A_prueba)>0.0001,2));
        a_store{jj} = A_prueba;
        jj = jj+1;
    end
    [tau_opt,ind] = min(costo);
    A_prueba = a_store{ind};

end

