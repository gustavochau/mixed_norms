function Z = zupdate_0inf( X,V,rho,lambda,W)
    A = W.*X+V;
    [Z, tau_opt,iter] = prox_newton_pruned( A,lambda/rho);
end