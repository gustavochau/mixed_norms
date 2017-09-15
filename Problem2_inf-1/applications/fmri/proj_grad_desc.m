function [ x, num_iter ,hist,tiempo] = proj_grad_desc( x0, proj_op, grad_op, tol, max_iter, alpha)
%PROJ_GRAD_DESC projected gradient descent routine. If step size alpha is not
%provided, it uses the spectal stepsizes of Barzilai and Borwein
%   Detailed explanation goes here
    if nargin<6
        bbflag = 1;
        alpha = 0.01;
    else
        bbflag = 0;
    end

    x_old = x0;
    x = x0;
    g_old = grad_op(x0);
    num_iter = max_iter;
    for ii=1:max_iter
        
        ii
        tic
        g = grad_op(x);
        
        t = x-alpha*g;
%         tic
        x = proj_op(t);
%         toc
        deltax=x(:)-x_old(:);
        hist(ii)=norm(deltax,2);
        if (norm(deltax,2)) < tol
            num_iter = ii;
            break
        end
%         if (bbflag && ii>1)
%             deltag=g(:)-g_old(:);
%             alpha = (deltax'*deltax)/(deltag'*deltax);
%             g_old = g;
%         end
        x_old = x;
        tiempo(ii)=toc;
    end
    

end

