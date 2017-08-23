clc;
clear;
close all;

k = 1;
f_t = @(u) k*(u.^2-1).*(u<1) + (3*k*(u.^2-1)).*(u>=1);
% f_t = @(u) k*(u.^4+u).*(u<0) -k* (u.^3+u).*(u>=0);
x = -5:0.1:5;
plot(x,f_t(x));

x0=-9;
tol = 1E-7;
max_iter = 4000;

% error = abs(f_t(x));
% iter = 1;
% while (error>tol)
%     x_hist(iter) = x;
% 
%     Fx = f_t(x);
%     error = abs(Fx);
%     if error < tol
%         disp('tol')
%         break;
%     end
%     Fy = f_t(x+Fx);
%     delta = (Fy-Fx)/(Fx);
%     x = x - Fx/delta;
%     
%     if (iter>= max_iter)    
%         break
%     end
%     iter = iter+1;
% end

% [x_hist,error]=steff_aitken(f_t,x0,tol,max_iter);
  i=0;
  err=tol+1;
  x=x0;
  phi=0;
  while(i<max_iter && err>tol)
      x_hist(i+1) = x;
      xx=x;
      fxk=f_t(x);
%       tolf=tol*abs(phi);
      if abs(fxk)<=tol
          disp('tol')
         break
      end
      fxk2=f_t(x+fxk);
      phi=(fxk2-fxk)/fxk;
      x=xx-fxk/phi;
      err=abs(x-xx);
      i=i+1;
  end
abs(x)
i

[x_a,x_a_hist,error_a,iter_a] = steff_amat(f_t,x0,tol,max_iter,1E-9);
iter_a
abs(x_a)
