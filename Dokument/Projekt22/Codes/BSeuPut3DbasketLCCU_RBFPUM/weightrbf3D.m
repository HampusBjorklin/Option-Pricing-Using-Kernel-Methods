function [phi,phix,phiy,phiz,phixx,phiyy,phizz,phixy,phixz,phiyz] = weightrbf3D(xc,locpts,xp,R)

Np = size(xp,1);
for i=1:length(Np)
    N(i) = length(locpts(i).ind);
end
zs = spalloc(size(xc,1),Np,sum(N));

[phi,phix,phiy,phiz,phixx,phiyy,phizz,phixy,phixz,phiyz] = deal(zs);
for i = 1:Np
    rho = 1/R;
    dx = xc(locpts(i).ind,1)-xp(i,1);
    dy = xc(locpts(i).ind,2)-xp(i,2);
    dz = xc(locpts(i).ind,3)-xp(i,3);
    r = sqrt(dx.^2 + dy.^2 + dz.^2);
    pos = find(r<=eps);
    r(pos)=r(pos) + eps; % To avoid division by zero

    phi(locpts(i).ind,i) = (4*rho*r+1).*(1-rho*r).^4;

    phix(locpts(i).ind,i) = -20*rho^2*dx.*(1-rho*r).^3;
    phiy(locpts(i).ind,i) = -20*rho^2*dy.*(1-rho*r).^3;
    phiz(locpts(i).ind,i) = -20*rho^2*dz.*(1-rho*r).^3;

    phixx(locpts(i).ind,i) = -20*rho*(1-rho*r).^2 + 60*rho^5*dx.^2.*r ...
        - 120*rho^4*dx.^2 + 60*rho^3*dx.^2./r;
    phiyy(locpts(i).ind,i) = -20*rho*(1-rho*r).^2 + 60*rho^5*dy.^2.*r ...
        - 120*rho^4*dy.^2 + 60*rho^3*dy.^2./r;
    phizz(locpts(i).ind,i) = -20*rho*(1-rho*r).^2 + 60*rho^5*dz.^2.*r ...
        - 120*rho^4*dz.^2 + 60*rho^3*dz.^2./r;

    phixy(locpts(i).ind,i) = 60*rho^5*dx.*dy.*r - 120*rho^4*dx.*dy ...
        + 60*rho^3*dx.*dy./r;
    phixz(locpts(i).ind,i) = 60*rho^5*dx.*dz.*r - 120*rho^4*dx.*dz ...
        + 60*rho^3*dx.*dz./r;
    phiyz(locpts(i).ind,i) = 60*rho^5*dy.*dz.*r - 120*rho^4*dy.*dz ...
        + 60*rho^3*dy.*dz./r;
end




