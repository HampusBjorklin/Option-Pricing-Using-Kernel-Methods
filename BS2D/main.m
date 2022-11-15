%
% Problem parameters
%
K=15; % Strike
T = 1; % Length of contract
r = 0.02; % Interest rate
sig1 = 0.2; % Volatility asset 1
sig2 = 0.15; % Volatility asset 2
rho = 0.5; % Correlation
N = 60; % Number of asset steps in one dimensions
M = 30; % Number of time steps
ep = 10; % Shape parameter (relate also to number of patches)
Np = 4; % 4x4 patches
del = 0.15; % Relative patch overlap

%
% We pick evaluation points on a grid for plotting
%
x = linspace(0,2.5*K,21);
[xx,yy] = meshgrid(x);
S = [xx(:), yy(:)];

U = BSeuCall2D_RBFPUM(S,K,T,r,sig1,sig2,rho,N,M,ep,Np,del);

%
% Plot the solution
%
uu = reshape(U,size(xx));

H = mesh(xx,yy,uu);
set(H,'FaceColor','interp')
axis equal
%
% Plot the difference between the solution and the payoff, although the payoff is now defined lower....
%
payoff = @(S,K,r,t) max(0.5*(S(:,1)+S(:,2))-K*exp(-r*t),0);
diffU = U - payoff(S,K,r,T);
diffuu = reshape(diffU,size(xx));

figure
H = mesh(xx,yy,diffuu);
set(H,'FaceColor','interp')
axis equal
