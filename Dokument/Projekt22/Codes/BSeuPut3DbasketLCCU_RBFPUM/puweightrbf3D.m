function pu = puweightrbf3D(xc,locpts,xp,R)

Np = size(xp,1);
[phi,phix,phiy,phiz,phixx,phiyy,phizz,phixy,phixz,phiyz] = weightrbf3D(xc,locpts,xp,R);


% Compute the sums of the generating functions and their derivatives
s = sum(phi,2);

sx = sum(phix,2);
sy = sum(phiy,2);
sz = sum(phiz,2);

sxx = sum(phixx,2);
syy = sum(phiyy,2);
szz = sum(phizz,2);

sxy = sum(phixy,2);
sxz = sum(phixz,2);
syz = sum(phiyz,2);

for i=1:Np
    loc = locpts(i).ind;
    if (length(loc>0))
    pu(i).w = phi(loc,i)./s(loc);

    pu(i).wx = phix(loc,i)./s(loc) - phi(loc,i).*sx(loc)./s(loc).^2;
    pu(i).wy = phiy(loc,i)./s(loc) - phi(loc,i).*sy(loc)./s(loc).^2;
    pu(i).wz = phiz(loc,i)./s(loc) - phi(loc,i).*sz(loc)./s(loc).^2;

    pu(i).wxx = -2*phix(loc,i).*sx(loc)./s(loc).^2 + phixx(loc,i)./s(loc) ...
        + phi(loc,i).*(2*sx(loc).^2./s(loc).^3 - sxx(loc)./s(loc).^2);
    pu(i).wyy = -2*phiy(loc,i).*sy(loc)./s(loc).^2 + phiyy(loc,i)./s(loc) ...
        + phi(loc,i).*(2*sy(loc).^2./s(loc).^3 - syy(loc)./s(loc).^2);
    pu(i).wzz = -2*phiz(loc,i).*sz(loc)./s(loc).^2 + phizz(loc,i)./s(loc) ...
        + phi(loc,i).*(2*sz(loc).^2./s(loc).^3 - szz(loc)./s(loc).^2);
    
    pu(i).wxy = phixy(loc,i)./s(loc) - phix(loc,i).*sy(loc)./s(loc).^2 ...
        - phiy(loc,i).*sx(loc)./s(loc).^2 - phi(loc,i).*sxy(loc)./s(loc).^2 ...
        + 2*phi(loc,i).*sx(loc).*sy(loc)./s(loc).^3;
    pu(i).wxz = phixz(loc,i)./s(loc) - phix(loc,i).*sz(loc)./s(loc).^2 ...
        - phiz(loc,i).*sx(loc)./s(loc).^2 - phi(loc,i).*sxz(loc)./s(loc).^2 ...
        + 2*phi(loc,i).*sx(loc).*sz(loc)./s(loc).^3;
    pu(i).wyz = phiyz(loc,i)./s(loc) - phiy(loc,i).*sz(loc)./s(loc).^2 ...
        - phiz(loc,i).*sy(loc)./s(loc).^2 - phi(loc,i).*syz(loc)./s(loc).^2 ...
        + 2*phi(loc,i).*sy(loc).*sz(loc)./s(loc).^3; 
    else
      pu(i).w=[];
    end
    
    end

  
  
  
  
  