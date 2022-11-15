%Copyright (C) 2015 Victor Shcherbakov
% Modified Elisabeth Larsson 2022

    %This file is part of BENCHOP.
    %BENCHOP is free software: you can redistribute it and/or modify
    %it under the terms of the GNU General Public License as published by
    %the Free Software Foundation, either version 3 of the License, or
    %(at your option) any later version.

    %BENCHOP is distributed in the hope that it will be useful,
    %but WITHOUT ANY WARRANTY; without even the implied warranty of
    %MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %GNU General Public License for more details.

    %You should have received a copy of the GNU General Public License
    %along with BENCHOP. If not, see <http://www.gnu.org/licenses/>.

function U=EuRBFPU2D(S,K,T,r,sig1,sig2,rho,payoff,nearbc,neard,farbc,fard,phi,ep,xc,M,npx,npy,delta)

% Compute the value of a European option using uniform RBF-PUM approximation
% S(1:Ne,1:2) evaluation points, points we want to know the price at
% K strike
% T maturity time
% r risk free interest rate
% sig1 volatility on 1st asset x direction
% sig2 volatility on 2nd asset y direction
% rho correlation
% payoff(S,K,r,t) function handle
% nearbc(S,K,r,t) function handle boundary condition at zero normally
% neard the distance below which we apply the near bc  
% farbc(S,K,r,t) function handle boundary condition far from strike
% fard the distance above which we apply the far bc    
% phi is the RBF to be used e.g. 'mq', 'gs'
% ep is the (constant) shape parameter  
% xc(1:N,1:2) are center points
% M is the number of time steps
% npx number of partitions in x direction
% npy number of partitions in y direction
% delta is the overlap parameter for the patches  

%  (C) Victor Shcherbakov & Elisabeth Larsson

% Define domain parameters
xeval = S; % evaluation points
Ne = size(S,1);

N = size(xc,1);

% Define boundary points
dx = sqrt(sum(xc.^2,2));
nearind = find(dx<=neard);
farind = find(dx>=fard);
inind = find(dx>neard & dx<fard);

% Define weights for BDF2 time stepping with fixed matrix
[k,beta0,beta1,beta2]=BDF2coeffs(T,M);

% Compute the size of the partitions and the patch centers
D = fard-neard;
Hx = D/npx;
x = neard+Hx/2:Hx:fard-Hx/2;
Hy = D/npy;
y = neard+Hy/2:Hy:fard-Hy/2;
[xx,yy] = meshgrid(x,y); % centers of partitions
xp = [xx(:) yy(:)];
Np = size(xp,1);
%
% Determine a patch radius such that one "box" is exactly enclosed
%
R0 = hypot(0.5*Hx,0.5*Hy);
R = R0*(1+delta);

[ptch(1:Np).R] = deal(R);
% Check which points go to which partition
for i=1:Np
    ptch(i).midpt = xp(i,:);
    flagin = hypot(xp(i,1)-xc(:,1),xp(i,2)-xc(:,2)) <= R;
    locpts(i).ind = find(flagin);
    Ni(i)=length(locpts(i).ind);
end

% Compute partition of unity weight functions at the center points
pu = puweightrbf(xc,locpts,ptch);

% Do the same for eval points
for i = 1:Np
    flagin = hypot(xp(i,1)-xeval(:,1),xp(i,2)-xeval(:,2)) <= R;
    loceval(i).ind = find(flagin);
    Nei(i)=length(loceval(i).ind);
end

% Compute partition of unity weight functions at the evaluation points
pu_e = puweightrbf(xeval,loceval,ptch);

derivativeflag = 1;
evalflag = 0;
% Form local RBF matrices and assemble the global B-S operator in L
L = spalloc(N,N,sum(Ni.^2));
for i=1:Np
    %
    % The operator is only computed for interior points
    %
    ind = locpts(i).ind;
    [~,locind] = intersect(locpts(i).ind,inind);
    % We always compute A0, but not necessarily the other parts
    xloc = xc(ind,:);
    rc = xcdist(xloc,xloc,derivativeflag);
    A0  = RBFmat(phi,ep,rc,'0');
    if (length(locind)>0) % If there are interior points
      nloc = length(locind);
      Ax  = RBFmat(phi,ep,rc,'1',1);
      Ay  = RBFmat(phi,ep,rc,'1',2);
      Axx = RBFmat(phi,ep,rc,'2',1);
      Axy = RBFmat(phi,ep,rc,'m2',[1,2]);
      Ayy = RBFmat(phi,ep,rc,'2',2);
      %
      % Diagonal weight matrices
      %
      Wi = spdiags(pu(i).w(locind),0,nloc,nloc);
      Wix = spdiags(pu(i).wx(locind),0,nloc,nloc);
      Wiy = spdiags(pu(i).wy(locind),0,nloc,nloc);
      Wixx = spdiags(pu(i).wxx(locind),0,nloc,nloc);
      Wixy = spdiags(pu(i).wxy(locind),0,nloc,nloc);
      Wiyy = spdiags(pu(i).wyy(locind),0,nloc,nloc);
      %
      % Diagonal coefficient matrices
      %
      X = spdiags(xloc(locind,1),0,nloc,nloc);
      Y = spdiags(xloc(locind,2),0,nloc,nloc);
      %
      % Local differentiation matrices, computed using product rule on W*A
      % 
      D0 = Wi*A0(locind,:); 
      Dx = Wi*Ax(locind,:) + Wix*A0(locind,:);
      Dy = Wi*Ay(locind,:) + Wiy*A0(locind,:);

      Dxx = Wi*Axx(locind,:) + 2*Wix*Ax(locind,:) + Wixx*A0(locind,:); 
      Dyy = Wi*Ayy(locind,:) + 2*Wiy*Ay(locind,:) + Wiyy*A0(locind,:);

      Dxy = Wi*Axy(locind,:) + Wix*Ay(locind,:) + ...
	    Wiy*Ax(locind,:) + Wixy*A0(locind,:);
      % 
      % Evaluate the Black-Scoles operator locally. Note that we multiply by
      % the inverse of A0 to change form coefficient to nodal representation.
      %
      B = (r*X*Dx + r*Y*Dy ...
           + 0.5*sig1^2*X.^2*Dxx ...
           + rho*sig1*sig2*X*Y*Dxy ...
           + 0.5*sig2^2*Y.^2*Dyy - r*D0)/A0;
      %
      % Add the contributions to the global matrix
      %
      L(ind(locind),locpts(i).ind) = L(ind(locind),locpts(i).ind) + B;
    end
    ind_e = loceval(i).ind;
    %
    % Form the local evaluation matrix (here we could also evaluate greeks)
    % 
    if length(ind_e)>0
        re = xcdist(xeval(ind_e,:),xloc,evalflag);
        Wei = spdiags(pu_e(i).w,0,Nei(i),Nei(i));
        E{i} = Wei*RBFmat(phi,ep,re,'0')/A0; % evaluation matrix
    end

end

% Time-stepping matrix I-beta_0*L in the interior and I on the boundary
C = speye(N) - beta0*L;
%
% Use LU factorization for speed (not pivoted here, to preserve sparsity(?))
%
[Lc,Uc] = lu(C);
%
% Prepare for first time step. Note that the first time rhs is just the
% initial condition.
%
tn = k(1);
u0 = payoff(xc,K,r,0);
rhs = u0;
tvec = cumsum(k); % The time after each step


% Time-stepping loop
for n = 1:M
    %
    % New solution value from linear system
    % 
    u = Uc\(Lc\rhs);
    % The final step is never done, but we still compute
    % a right hand side
    nextstep = min(M,n+1); 
    
    rhs = beta1(nextstep)*u-beta2(nextstep)*u0;
    
    % boundary conditions at the next time level
    tn = tvec(nextstep); % Should be one step ahead        

    rhs(nearind) = nearbc(xc(nearind,:),K,r,tn);
    rhs(farind) = farbc(xc(farind,:),K,r,tn);
    
    % move the solutions one step
    u0 = u;
end

% Compute the solution or some derivative at the evaluation points
U = zeros(Ne,1);
for i=1:Np
  ind_e = loceval(i).ind;
  if (length(ind_e)>0)
    ind=locpts(i).ind;
    se=size(ind_e)
    si=size(ind)
    sE=size(E{i})
    U(ind_e) = U(ind_e) + E{i}*u(ind);
  end
end  










