function [U]=BSeuPut3DbasketLCCU_RBFPUM(tol)

T = 1; % Time to maturity
rho = [1, 0.25, 0.25;
       0.25, 1, 0.25;
       0.25, 0.25, 1];
sig = [0.2, 0.2, 0.2]; % Correlation matrix
r = 0.06; % Interest rate
K = 40; % Strike price
S = [40, 40, 40]; % Evaluation point
payoff = @(S,K,r,t) max(K*exp(-r*t)-sum(S,2)/3,0); % Payoff function


% Generate suitable nodes for the type of option and domain
ntype = 'uni'; % 'cluster' - clustered; 'uni' - uniform grid
mode = 'sinh'; % Type of clustering: 'sinh' or 'cheb'
ell = 2; % Amount of clustering
scale = K; % Scaling factor
focus = 3*K; % Focus of clustering

delta = 0.2; % Overlap size
dim = 3; % Problem dimensionality


if tol == 0.001 || tol == 10^(-3) || tol == 1e-3
    phi = 'mq'; % RBF type
    ep = 1.05; % Shape parameter
    N = 10; % # Space points
    M = 2; % # Time steps
    np = 3; % # Partitions
    % Generate suitable nodes for the type of option and domain
    xmin = 0; % Lower domain boundary
    xmax = 4*K; % Upper domain boundary   
elseif tol == 0.0001 || tol == 10^(-4) || tol == 1e-4
    phi = 'mq'; % RBF type
    ep = 1.05; % Shape parameter
    N = 14; % # Space points
    M = 7; % # Time steps
    np = 5; % # Partitions
    % Generate suitable nodes for the type of option and domain
    xmin = 0; % Lower domain boundary
    xmax = 4*K; % Upper domain boundary 
elseif tol == 0.00001 || tol == 10^(-5)  || tol == 1e-5
    phi = 'mq'; % RBF type
    ep = 1.01; % Shape parameter
    N = 14; % # Space points
    M = 15; % # Time steps
    np = 5; % # Partitions
    % Generate suitable nodes for the type of option and domain
    xmin = 0; % Lower domain boundary
    xmax = 4*K; % Upper domain boundary 
elseif tol == 0.000001 || tol == 10^(-6) || tol == 1e-6
    phi = 'mq'; % RBF type
    ep = 1; % Shape parameter
    N = 20; % # Space points
    M = 18; % # Time steps
    np = 5; % # Partitions
    % Generate suitable nodes for the type of option and domain
    xmin = 0; % Lower domain boundary
    xmax = 4*K; % Upper domain boundary 
end
    
xc = GetNodesSimplex(xmin,xmax,ntype,N,focus,ell,mode,dim);

U = BSeuPut3Dpum(S,K,T,r,sig,rho,payoff,phi,ep,xc,M,np,delta);

%U = U';
% Uref = 1.22309;
% err = abs(U-Uref);















