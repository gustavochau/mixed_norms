function [x,x_hist,error,iter]=root_newt(f,x0,tol,max_iter)

x=x0;
iter = max_iter;
for ii=1:max_iter
    x_hist(ii) = x;
  
    [Fx,DFx]=feval(f,x);
    error = abs(Fx);
    if ii>1
%         error = abs(x-x_hist(ii-1));
        if error < tol
            iter = ii;
            break
        end
    end
    
    x = x-Fx/DFx;
end

end