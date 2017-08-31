function [x,x_hist,error,iter]=root_newt(f,x0,tol,max_iter)

x=x0;
iter = max_iter;
for ii=1:max_iter
    x_hist(ii) = x;
  
    [Fx,DFx]=feval(f,x);
    error = abs(Fx);
    if error < tol
        iter = ii;
        break
    end
    x = x-Fx/DFx;
end

end