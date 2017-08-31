function [ f_val ] = f_prueba( tau,B,lambda )
%F_PRUEBA Summary of this function goes here
%   Detailed explanation goes here
    [nrows,ncols] = size(B);
    nb= sum(abs(B),2);
    for ii=1:nrows
        [x, loops(ii),gamma(ii)] = projL1Mich(B(ii,:)', tau, 20);
        A(ii,:) = x';
    end
    f_val=0;
    for ii=1:nrows
        if nb(ii) >= tau
%             if all(abs(B(ii,:))<gamma(ii))
%                 f_val = f_val + max(abs(B(ii,:)));
%             else
%                 f_val = f_val + gamma(ii);
%             end

%             f_val= f_val + max(abs(B(ii,:))-max(abs(B(ii,:))-gamma(ii),0));
            
              z =  sign(B(ii,:)').*double(abs(B(ii,:))'>gamma(ii));
              f_val = f_val + ((z')*(B(ii,:)')-tau)/((z')*z);
        end
    end
    f_val = f_val-lambda;
end

