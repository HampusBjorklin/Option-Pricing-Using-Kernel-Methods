%
%  2D Solver of Eu call option using Holgers method of reproducing Kernel
%   Solved in rotated system! Lite fulkod men it works okay!
% Something i off close to the boundary. Seems to be with the PDE Solver
function [U, u, X, XT] = Holger2DEuCallTransform(X_eval, smax, K, T, r, sig1, sig2, rho, anchor, n, M, ep)
%% Evaluation and Center Points
% Choose a grid Omega in [1/3 K, 5/3*K]^d to evalutate
dim = 2;
Ne = length(X_eval);    %Number of evaluation points
N = dim*n + 1; % Total number of centerpoints


anchor = anchor/smax;
% anchor = [0.5, 0.5]
% %Oscar cool way to fix the center points.
domain = [0,1];
x = linspace(0,1, n);


% Scale down the Evalutaion points to region [0, 1]
X_eval = X_eval./smax;
K = K /smax;

% Mat = GCV2S(dim);
%Transform functions
f_s2v = @(S) [(S(:,1) + S(:,2))/2, ...
    (S(:,1) - S(:,2))/2];

f_v2s = @(V)  [V(:,1) + V(:,2), ...
    V(:,1) - V(:,2)];



% f_v2s = @(V) (Mat*V')';
% f_s2v = @(S) (Mat'*S')';
% f_s2v = @(X) [(X(:,1) + X(:,2))/2, ...
%                 (X(:,1) - X(:,2) + 1)/2];

X = getXVector(f_s2v(anchor'), n);
XT = f_v2s(X);
X_eval = f_s2v(X_eval);
% Define boundary points
distClose = 0;
distFar = 1;
range = sqrt(sum(XT.^2,2));
indClose = find(range<=distClose);
indFar = find(range>=distFar);
indInter = find(range>distClose & range<distFar);

%Obtain the set of center points

% Boundary Conditions
closeBC = @(S,K,r,t) zeros(size(S,1),1);
farBC   = @(S,K,r,t) max(0.5*(S(:,1)+S(:,2))-K*exp(-r*t),0);




%Find interior points
XInter = X(indInter,:);
NInter = length(indInter);


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

B0 = A0(indInter, :);
Bx = 1/2.*(Ax(indInter, :) + Ay(indInter, :));
By = 1/2.*(Ax(indInter, :) - Ay(indInter, :));
Bxx = 1/2.*(1/2.*Axx(indInter, :) + Axy(indInter, :) + 1/2.*Ayy(indInter, :));
Bxy = 1/4*(Axx(indInter, :) - Ayy(indInter, :));
Byy = 1/2.*(1/2.*Axx(indInter, :) - Axy(indInter, :) + 1/2.*Ayy(indInter, :));


s1 = XT(:,1);
s2 = XT(:,2);
S1 = spdiags(s1(indInter,1),0,NInter,NInter);
S2 = spdiags(s2(indInter,1),0,NInter,NInter);



Operator = (r*S1*Bx + r*S2*By ...
    + 0.5*sig1^2*S1.^2*Bxx ...
    + rho*sig1*sig2*S1*S2*Bxy ...
    + 0.5*sig2^2*S2.^2*Byy - r*B0)/A0;
% Operator should be R^N_inner x N
%Add Interiour operator at correct point in the Large operator
L = zeros(N,N);
L(indInter, :) = L(indInter, :) + Operator;

%% Solver
% Using the BDF2 function from Elisabeths kod.

k = T/(M-1);
beta0 = k* 2/3;
beta1 = 4/3;
beta2 = 1/3;

C = eye(N) -beta0*L;

%Initialization
u0 = farBC(XT, K, r, 0); %IC
rhs = u0;
tvec = 0:k:T;

[Lc,Uc] = lu(C); %Tips?

%Time steppin'
if T == 0
    tvec = zeros(M,1);
end
for m = 1:M
    u =Uc\(Lc\rhs);
    
    % BC at the next time level
    nextstep = min(M,m+1);
    rhs = beta1*u - beta2*u0;
    
    tn = tvec(nextstep); % Should be one step ahead
    
    rhs(indClose,:) = closeBC(XT(indClose,:),K,r,tn);
    rhs(indFar,:) = farBC(XT(indFar,:),K,r,tn);
    
    % move the solutions one step
    u0 = u;
end


%% Obtain solution at Desired Evalutation Points
% Vet ej om randpunkterna ska vara med n�r vi best�mmer eval.
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