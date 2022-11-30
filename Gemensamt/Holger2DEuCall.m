%
%  2D Solver of Eu call option using Holgers method of reproducing Kernel
%
function [U, u, X] = Holger2DEuCall(X_eval, smax, K, T, r, sig1, sig2, rho, anchor, n, M, ep)
%% Evaluation and Center Points
% Choose a grid Omega in [1/3 K, 5/3*K]^d to evalutate
dim = 2;
Ne = length(X_eval);    %Number of evaluation points
N = dim*n + 1; % Total number of centerpoints 


% Scale down the Evalutaion points to region [0, 1]
X_eval = X_eval./smax;
K = K /smax;
anchor = anchor/smax;

%Obtain the set of center points
% X = [[x',anchor*ones(n,1)]; [anchor, anchor]; [anchor*ones(n,1), x']];
X = getXVector(anchor, n, 0);

% Define boundary points
distClose = 0;
distFar = 1;
range = sqrt(sum(X.^2,2));
indClose = find(range<=distClose);
indFar = find(range>=distFar);
indInter = find(range>distClose & range<distFar);

% Boundary Conditions
closeBC = @(S,K,r,t) zeros(size(S,1),1);
farBC   = @(S,K,r,t) max(0.5*(S(:,1)+S(:,2))-K*exp(-r*t),0);
   


%Find interior points
XInter = X(indInter,:);
NInter = length(indInter);
XX = spdiags(X(:,1),0,N,N);
YY = spdiags(X(:,2),0,N,N);


%Local Differentiation matrices
A0 = zeros(N, N); Ax = zeros(N, N);
Ay = zeros(N, N); Axx = zeros(N, N);
Axy = zeros(N, N); Ayy = zeros(N, N);

eps2 = ep^2;
for i = 1:N
    for j = 1:N
        A0(i, j) = RepKernel(X(i,:), X(j, :), eps2);
        Ax(i, j) = Dx1RepKernel(X(i,:), X(j,:), eps2);
        Ay(i, j) = Dx2RepKernel(X(i,:), X(j,:), eps2);
        Axx(i, j) = Dxx1RepKernel(X(i,:), X(j,:), eps2);
        Axy(i, j) = Dxx12RepKernel(X(i,:), X(j,:), eps2);  %This is zero in this case eksde
        Ayy(i, j) = Dxx2RepKernel(X(i,:), X(j,:), eps2);
    end
end

%Black schouls(?) Operator (This is for interiour points.
Operator = (r*XX*Ax + r*YY*Ay ...
           + 0.5*sig1^2*XX.^2*Axx ...
           + rho*sig1*sig2*XX*YY*Axy ...
           + 0.5*sig2^2*YY.^2*Ayy - r*A0)/A0;
       
%Add Interiour operator at correct point in the Large operator
L = spalloc(N,N,sum(NInter.^2));
% L(indInter, indInter) = L(indInter, indInter) + Operator;
L = Operator;

%% Solver
% Using the BDF2 function from Elisabeths kod. 
[k,beta0,beta1,beta2]=BDF2coeffs(T,M);
C = speye(N) - beta0*L; 
 
%Initialization
u0 = farBC(X, K, r, 0); %IC
tn = k(1);
rhs = u0;
tvec = cumsum(k);

%Time steppin'
for m = 1:M
    u =C\rhs;
    
    % BC at the next time level
    nextstep = min(M,m+1);
    tn = tvec(nextstep); % Should be one step ahead
    rhs = beta1(nextstep)*u - beta2(nextstep)*u0;
    rhs(indClose,:) = closeBC(X(indClose,:),K,r,tn);
    rhs(indFar,:) = farBC(X(indFar,:),K,r,tn);
    
    % move the solutions one step
    u0 = u;
end


%% Obtain solution at Desired Evalutation Points
% Vet ej om randpunkterna ska vara med när vi bestämmer eval.
% Evaluation matrix
E = zeros(Ne, N);
for i = 1:Ne
    for j = 1:N
        E(i, j) = RepKernel(X_eval(i,:), X(j,:), eps2);
    end
end
% E = E/A0; %Want to transform from nodal representation
E = E/A0;
% Solution
U = E*u;
%Rescale 
X_eval = X_eval*smax;
U = U*smax; K = K*smax;
end

%% All the kernel funcitons
function k = RepKernel(x, y, eps2)
    k = 1 + sqrt(1 + eps2*(x(1) - y(1))^2) + sqrt(1 + eps2*(x(2) - y(2))^2);
end
function k = Dx1RepKernel(x, y, eps2)
    
    k = eps2*(x(1) - y(1))/ (sqrt(eps2*(x(1) - y(1))^2 + 1));
end
function k = Dx2RepKernel(x, y, eps2)
    
    k = eps2*(x(2) - y(2))/ (sqrt(eps2*(x(2) - y(2))^2 + 1));
end
function k = Dxx1RepKernel(x, y, eps2)
    
    k = eps2/ ((eps2*(x(1) - y(1))^2 + 1)^(3/2));
end
function k = Dxx2RepKernel(x, y, eps2)
    
    k = eps2/ ((eps2*(x(2) - y(2))^2 + 1)^(3/2));
end
function k = Dxx12RepKernel(x, y, eps2)
    k = 0;
end