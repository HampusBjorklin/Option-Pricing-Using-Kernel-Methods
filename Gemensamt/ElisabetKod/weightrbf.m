% Evaluate Wendland's C² function and its derivatives

%  (C) Alfa Heryudono 2010, Victor Shcherbakov, Elisabeth Larsson 2015

function [phi,phix,phiy,phixx,phixy,phiyx,phiyy] = weightrbf(xc,locpts,xp,R)

  Np = size(xp,1);
  for i=1:length(Np)
    N(i) = length(locpts(i).ind);
  end
  zs = spalloc(size(xc,1),Np,sum(N));
  
  [phi,phix,phiy,phixx,phixy,phiyx,phiyy] = deal(zs);
  for i = 1:Np
    rho = 1/R(i);
    dx = xc(locpts(i).ind,1)-xp(i,1);
    dy = xc(locpts(i).ind,2)-xp(i,2);
    r = sqrt(dx.^2 + dy.^2);
    pos = find(r<=eps);
    r(pos)=r(pos) + eps; % To avoid division by zero

    phi(locpts(i).ind,i) = (4*rho*r+1).*(1-rho*r).^4;
    
    phix(locpts(i).ind,i) = -20*rho^2*dx.*(1-rho*r).^3;
    phiy(locpts(i).ind,i) = -20*rho^2*dy.*(1-rho*r).^3;
    
    phixx(locpts(i).ind,i) = -20*rho*(1-rho*r).^2 + 60*rho^5*dx.^2.*r ...
        - 120*rho^4*dx.^2 + 60*rho^3*dx.^2./r;
    phixy(locpts(i).ind,i) = 60*rho^5*dx.*dy.*r - 120*rho^4*dx.*dy ...
        + 60*rho^3*dx.*dy./r;
    phiyx(locpts(i).ind,i) = 60*rho^5*dy.*dx.*r - 120*rho^4*dy.*dx ...
        + 60*rho^3*dy.*dx./r;
    phiyy(locpts(i).ind,i) = -20*rho*(1-rho*r).^2 + 60*rho^5*dy.^2.*r ...
        - 120*rho^4*dy.^2 + 60*rho^3*dy.^2./r;
end














