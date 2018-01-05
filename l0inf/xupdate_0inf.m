function X = xupdate_0inf( X,Z,V,B,rho,W)
    X = (rho*W.*(Z-V)+B)./(1+rho*W.*W);
end