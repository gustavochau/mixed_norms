function [x,x_hist,error,iter]=steff_amat_acc(f,x0,tol,max_iter,tol_u)

tol_c = 1E-12;
x=x0;
iter = max_iter;
for ii=1:max_iter
    x_hist(ii) = x;
  
    Fx=feval(f,x);
    error = abs(Fx);
    if error < tol
        iter = ii;
        break
    end
    alpha = mean(([tol_c tol_u]/(Fx*abs(Fx))));
    x2 = x + alpha*abs(Fx);
    Fy=feval(f,x2);
    delta = (Fy-Fx)/(alpha*abs(Fx));
    x = x-Fx/sqrt(delta);
end

end