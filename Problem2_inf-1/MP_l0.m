function [ x ] = MP_l0(y, tau )
%MP_L0 Summary of this function goes here
%   Detailed explanation goes here
    normas = y.^2;
    x = 0*y;
    [ordered,I] = sort(normas,'descend' );
    ind_sel = I(1:min(floor(tau),length(y)));
    x(ind_sel) = y(ind_sel);
end

