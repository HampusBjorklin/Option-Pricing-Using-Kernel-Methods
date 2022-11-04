function U = BSeuPut3Dpum(S,K,T,r,sig,rho,payoff,phi,ep,xc,M,np,delta)

% Scale the strike to reduce to unit-strike problem
K0 = K;
K = 1;

% Scale center points
xc = xc/K0;
Nc = size(xc,1);
dim = size(xc,2);

% Scale evaluation points
xe = S/K0;
Ne = size(xe,1);

% Find boundaries
xmin = min(xc(:,1));
xmax = max(xc(:,1));
bnear = find(xc(:,1)==xmin & xc(:,2)==xmin & xc(:,3)==xmin);
bfar = find(xc(:,1)+xc(:,2)+xc(:,3)>=xmax-0.01 & xc(:,1)+xc(:,2)+xc(:,3)<=xmax+0.01);
indb = [bnear; bfar];

% Find partition centers
[xm,Np,Rad] = GetPartitionCentersSimplex(xmin,xmax,np,delta,dim);

% Check which center points go to which partition
[ibox,Ni] = Point2PartitionBelonging(xc,xm,Rad,Np);

% Compute weight functions for each partition using Shepard's method
puw = puweightrbf3D(xc,ibox,xm,Rad);

% Check which evaluation points go to which partition
[ebox,Nei] = Point2PartitionBelonging(xe,xm,Rad,Np);

% Compute weight functions for each partition using Shepard's method
puw_e = puweightrbf3D(xe,ebox,xm,Rad);

% Allocate memory for the discrete operator matrix 
L = spalloc(Nc,Nc,sum(Ni.*Ni));

% Allocate memory for the evaluation matrix 
Eval = spalloc(Ne,Nc,sum(Ni.*Nei));

for i = 1:Np    
    % Find indeces 
    ind = ibox(i).ind;
    xloc = xc(ind,:);
       
    % Compute local RBF matrices 
    [B,A0] = BSop3D(xloc,puw(i),r,sig,rho,phi,ep);
    
    % Assemble 
    L(ind,ind) = L(ind,ind) + B;
     
    % Evaluation matrix
    ind_e = ebox(i).ind;
    flagempt = isempty(ind_e);
    if flagempt == 0
        re = xcdist(xe(ind_e,:),xloc,'0');
        Wei = spdiags(puw_e(i).w,0,Nei(i),Nei(i));
        E{i} = Wei*RBFmat(phi,ep,re,'0')/A0; % evaluation matrix
        Eval(ind_e,ind) = Eval(ind_e,ind) + E{i};
    end     
end
L(indb,:) = 0;

% Sparse identity
I = speye(Nc);

% BDF2 coefficients
[k,beta0,beta1,beta2]=BDF2coeffs(T,M);

% Black-Scholes operator
C = I - beta0*L;

% GMRES setup
tol = 1e-05;
maxit = 15;
setup.type = 'nofill';
[iL,iU] = ilu(C,setup);

% Initial condition
u0 = payoff(xc,K,r,0);
rhs = u0; 
tn = k(1);
rhs(indb) = payoff(xc(indb,:),K,r,tn);


% Time-stepping loop
for n = 1:M
    % Solve with GMRES
    [u,flag,relres,iter,resvec] = gmres(C,rhs,[],tol,maxit,iL,iU,u0);
    
    % Prepare the right hand side for the next step
    nextstep = min(M,n+1);
    tn = sum(k(1:nextstep)); % Should be one step ahead
    rhs = beta1(nextstep)*u-beta2(nextstep)*u0;
    
    % Boundary conditions
    rhs(indb) = payoff(xc(indb,:),K,r,tn);
    
    % Update
    u0 = u;   
end

U = Eval*u;
U = K0*U;
















