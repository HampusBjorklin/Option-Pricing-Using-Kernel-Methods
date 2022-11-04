% Generate the spatial part of the 3D Black-Scholes operator 
% for the partition of unity setup given local computational nodes, pu weights, 
% interest rate, volatility, shape parameter and type of basis function

function [B,A0] = BSop3D(x,puw,r,sig,rho,phi,ep)

    n = size(x,1);

    rc = xcdist(x,x,1);

    A0 = RBFmat(phi,ep,rc,'0');
    
    Ax = RBFmat(phi,ep,rc,'1',1);
    Ay = RBFmat(phi,ep,rc,'1',2);
    Az = RBFmat(phi,ep,rc,'1',3);
    
    Axx = RBFmat(phi,ep,rc,'2',1);
    Ayy = RBFmat(phi,ep,rc,'2',2);
    Azz = RBFmat(phi,ep,rc,'2',3);
    
    Axy = RBFmat(phi,ep,rc,'m2',[1,2]);
    Axz = RBFmat(phi,ep,rc,'m2',[1,3]);
    Ayz = RBFmat(phi,ep,rc,'m2',[2,3]);

    % Compute matrices with PU weight on diagonal
    W = spdiags(puw.w,0,n,n);
    
    Wx = spdiags(puw.wx,0,n,n);
    Wy = spdiags(puw.wy,0,n,n);
    Wz = spdiags(puw.wz,0,n,n);
    
    Wxx = spdiags(puw.wxx,0,n,n);
    Wyy = spdiags(puw.wyy,0,n,n);
    Wzz = spdiags(puw.wzz,0,n,n);
    
    Wxy = spdiags(puw.wxy,0,n,n);
    Wxz = spdiags(puw.wxz,0,n,n);
    Wyz = spdiags(puw.wyz,0,n,n);
    
    % Compute diagonal matrices for the valriables
    X = spdiags(x(:,1),0,n,n);
    Y = spdiags(x(:,2),0,n,n);
    Z = spdiags(x(:,3),0,n,n);

    % Compute RBF-PUM differentiation matrices
    D0 = W*A0;

    Dx = W*Ax + Wx*A0;
    Dy = W*Ay + Wy*A0;
    Dz = W*Az + Wz*A0;

    Dxx = W*Axx + 2*Wx*Ax + Wxx*A0;
    Dyy = W*Ayy + 2*Wy*Ay + Wyy*A0;
    Dzz = W*Azz + 2*Wz*Az + Wzz*A0;

    Dxy = W*Axy + Wx*Ay + Wy*Ax + Wxy*A0;
    Dxz = W*Axz + Wx*Az + Wz*Ax + Wxz*A0;
    Dyz = W*Ayz + Wy*Az + Wz*Ay + Wyz*A0;

    B = (r*X*Dx ...
        + r*Y*Dy ...
        + r*Z*Dz ...
        - r*D0 ...
        + 0.5*sig(1)^2*X.^2*Dxx ...
        + 0.5*sig(2)^2*Y.^2*Dyy ...
        + 0.5*sig(3)^2*Z.^2*Dzz ...
        + rho(1,2)*sig(1)*sig(2)*X*Y*Dxy ...
        + rho(1,3)*sig(1)*sig(3)*X*Z*Dxz ...
        + rho(2,3)*sig(2)*sig(3)*Y*Z*Dyz)/A0;









