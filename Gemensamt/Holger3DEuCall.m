%
%  2D Solver of Eu call option using Holgers method of reproducing Kernel
%
function [U, u, X] = Holger3DEuCall(X_eval, smax, K, T, r, C, anchor, n, M, ep)
%% Evaluation and Center Points
% Choose a grid Omega in [1/3 K, 5/3*K]^d to evalutate
dim = 3;
Ne = length(X_eval);    %Number of evaluation points
N = dim*n + 1; % Total number of centerpoints 

% Scale down the Evalutaion points to region [0, 1]
X_eval = X_eval./smax;
K = K /smax;
anchor = anchor/smax;


%Transform functions
% Rotate = @(S) [(S(:,1) + S(:,2))/2, (1 + S(:,1) - S(:,2))/2];
X = getXVector(anchor, n);

% f_s2v = @(S) [(S(:,1) + S(:,2))/2, ...
%                 (S(:,1) - S(:,2))/2];
%             
% f_v2s = @(V)  [V(:,1) + V(:,2), ...
%                 V(:,1) - V(:,2)];

f_v2s = @(V) ([0.5, 1, 0.5; ...
              -0.5, 1, -0.5;...
              -1, 0, 1] * (V'))';
f_s2v = @(S) 1/2.*([1, -1, -1;...
                    1, 1, 0;...
                    1, -1, 1] * (S'))';

XT = fI(X);
X_eval = f(X_eval);


% Define boundary points
distClose = 0;
distFar = 1;
range = sqrt(sum(X.^2,2));
indClose = find(range<=distClose);
indFar = find(range>=distFar);
indInter = find(range>distClose & range<distFar);

% Boundary Conditions
closeBC = @(S,K,r,t) zeros(size(S,1),1);
farBC   = @(S,K,r,t) max(1/3*(S(:,1)+S(:,2)+S(:,3))-K*exp(-r*t),0);
   


%Find interior points
XInter = X(indInter,:);
NInter = length(indInter);
XX = spdiags(X(indInter,1),0,NInter,NInter);
YY = spdiags(X(indInter,2),0,NInter,NInter);
ZZ = spdiags(X(indInter,2),0,NInter,NInter);


%Local Differentiation matrices
A0 = zeros(N, N); Ax = zeros(N, N);
Ay = zeros(N, N); Axx = zeros(N, N);
Axy = zeros(N, N); Ayy = zeros(N, N);
eps2 = ep^2;
for i = 1:NInter
    for j = 1:NInter
        A0(i, j) = RepKernel(XInter(i,:), XInter(j, :), eps2);
        Ax(i, j) = DxRepKernel(XInter(i,:), XInter(j,:), eps2, 1);
        Ay(i, j) = DxRepKernel(XInter(i,:), XInter(j,:), eps2, 2);
        Az(i, j) = DxRepKernel(XInter(i,:), XInter(j,:), eps2, 3);
        Axx(i, j) = DxxRepKernel(XInter(i,:), XInter(j,:), eps2, 1);
        Ayy(i, j) = DxxRepKernel(XInter(i,:), XInter(j,:), eps2, 2);
        Azz(i, j) = DxxRepKernel(XInter(i,:), XInter(j,:), eps2, 3);
        Axy(i, j) = DxyRepKernel(XInter(i,:), XInter(j,:), eps2, 1, 2);
        Axz(i, j) = DxyRepKernel(XInter(i,:), XInter(j,:), eps2, 1, 3);
        Ayz(i, j) = DxyRepKernel(XInter(i,:), XInter(j,:), eps2, 2, 3);
    end
end

%Black schouls(?) Operator (This is for interiour points.
Operator = (r*XX*Ax + r*YY*Ay + r*ZZ*Az...
           + 0.5*C(1,1)*XX.^2*Axx ...
           + 0.5*C(2,2)*YY.^2*Ayy ...
           + 0.5*C(3,3)*YY.^2*Azz ...
           + 0.5*(C(1,2) +  C(2,1))*XX*YY*Axy ...
           + 0.5*(C(1,3) +  C(3,1))*XX*ZZ*Axz ...
           + 0.5*(C(2,3) +  C(3,2))*YY*ZZ*Axy ...
           - r*A0)/A0;
       
%Add Interiour operator at correct point in the Large operator
L = spalloc(N,N,sum(NInter.^2));
L(indInter, indInter) = L(indInter, indInter) + Operator;

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
E = zeros(Ne, NInter);
% Atot = zeros(N,N);
% for i = 1:N
%     for j = 1:N
%         Atot(i,j) = RepKernel(X(i,:), X(j,:), eps2);
%     end
% end
for i = 1:Ne
    for j = 1:NInter
        E(i, j) = RepKernel(X_eval(i,:), XInter(j,:), eps2);
    end
end
E = E/A0; %Want to transform from nodal representation

% Solution
U = E*u(indInter);
%Rescale 
X_eval = X_eval*smax;
U = U*smax; K = K*smax;
end

%% All the kernel funcitons
function k = RepKernel(x, y, eps2)
    k = 1 + sqrt(1 + eps2*(x(1) - y(1))^2) ...  % 1 + k1
        + sqrt(1 + eps2*(x(2) - y(2))^2)   ... % k2
        + sqrt(1 + eps2*(x(3) - y(3))^2)   ... % k3
        + sqrt(1 + eps2*(x(1) - y(1))^2)*sqrt(1 + eps2*(x(2) - y(2))^2) ... %k1*k2
        + sqrt(1 + eps2*(x(1) - y(1))^2)*sqrt(1 + eps2*(x(3) - y(3))^2) ... %k1*k3
        + sqrt(1 + eps2*(x(2) - y(2))^2)*sqrt(1 + eps2*(x(3) - y(3))^2);    %k2*k3
end
function k = DxRepKernel(x, y, eps2, dim)
    index = [1:3];
    index(index == dim) = [];
    
    k = eps2*(x(dim) - y(dim))/ (sqrt(eps2*(x(dim) - y(dim))^2 + 1)) ...
        *( 1 + sqrt(1 + eps2*(x(index(1)) - y(index(1)))^2) + sqrt(1 + eps2*( x(index(2)) - y(index(2)) )^2));
    
end
function k = DxxRepKernel(x, y, eps2, dim)
    index = [1:3];
    index(index == dim) = [];
    
    k = eps2/ ((eps2*(x(dim) - y(dim))^2 + 1)^(3/2)) ...
        *( 1 + sqrt(1 + eps2*(x(index(1)) - y(index(1)))^2) + sqrt(1 + eps2*(x(index(2)) - y(index(2)))^2));
end
function k = DxyRepKernel(x, y, eps2, dim1, dim2)
    k = eps2^2 * ((x(dim1) - y(dim1))*(x(dim2) - y(dim2))) ...
    /(sqrt(1 + eps2 *(x(dim1) - y(dim1))^2) * sqrt(1 + eps2 *(x(dim2) - y(dim2))^2));
end