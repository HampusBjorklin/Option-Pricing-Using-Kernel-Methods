% Evaluate partition of unity weight functions using Shepard's method
%
% (C) Elisabeth Larsson 2015
% Based on code initially developed by Alfa Heryudono, modified by Victor
% Shcherbakov and Elisabeth Larsson

function pu = puweightrbf(xc,locpts,ptch)
Np = length(ptch);
for k=1:Np
  R(k)=ptch(k).R;
  xp(k,:)=ptch(k).midpt;
end

[phi,phix,phiy,phixx,phixy,phiyx,phiyy] = weightrbf(xc,locpts,xp,R);
%
% Compute the sums of the generating functions and their derivatives
%
s = sum(phi,2);

sx = sum(phix,2);
sy = sum(phiy,2);

sxx = sum(phixx,2);
sxy = sum(phixy,2);
syx = sum(phiyx,2);
syy = sum(phiyy,2);

for i=1:Np
  loc = locpts(i).ind;
  pu(i).w = phi(loc,i)./s(loc);
  
  pu(i).wx = phix(loc,i)./s(loc) - phi(loc,i).*sx(loc)./s(loc).^2;
  pu(i).wy = phiy(loc,i)./s(loc) - phi(loc,i).*sy(loc)./s(loc).^2;

  pu(i).wxx = -2*phix(loc,i).*sx(loc)./s(loc).^2 + phixx(loc,i)./s(loc) ...
      + phi(loc,i).*(2*sx(loc).^2./s(loc).^3 - sxx(loc)./s(loc).^2);
  pu(i).wxy = phixy(loc,i)./s(loc) - phix(loc,i).*sy(loc)./s(loc).^2 ...
      - phiy(loc,i).*sx(loc)./s(loc).^2 - phi(loc,i).*sxy(loc)./s(loc).^2 ...
      + 2*phi(loc,i).*sx(loc).*sy(loc)./s(loc).^3;
  pu(i).wyx = phiyx(loc,i)./s(loc) - phiy(loc,i).*sx(loc)./s(loc).^2 ...
      - phix(loc,i).*sy(loc)./s(loc).^2 - phi(loc,i).*syx(loc)./s(loc).^2 ...
      + 2*phi(loc,i).*sy(loc).*sx(loc)./s(loc).^3; 
  pu(i).wyy = -2*phiy(loc,i).*sy(loc)./s(loc).^2 + phiyy(loc,i)./s(loc) ...
      + phi(loc,i).*(2*sy(loc).^2./s(loc).^3 - syy(loc)./s(loc).^2);

end



