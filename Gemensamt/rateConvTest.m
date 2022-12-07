

%% 2D
close all; clear all; clc;
%Ecomonic Parameters
K=  20; % Strike
T = 1; % Length of contract
r = 0.02; % Interest rate
sig1 = 0.2; % Volatility asset 1
sig2 = 0.1; % Volatility asset 2
rho = 0.5; % Correlation


%Numerical Parameters
dim = 2;
n_vector = 5:25; %Points in each dimention
M = 50; %Number of timesteps.



% Points for evaluation 
smax = 4*K;        %Largets value for simulation (center points)
Eval_smin = 0;  %Evalutaion min
Eval_smax = 4*K;  %Evaluation max
temp_x = linspace(Eval_smin,Eval_smax,41);
[xx, yy] = meshgrid(temp_x);
X_eval = [xx(:) yy(:)]; %Evaluation points 

ep = 15; % Shape parameter
anchor = [40; 40]; % Anchor, freezing point
% tic
% [U,u, X] = Holger2DEuCall(X_eval,smax, K, T, r, sig1,sig2, rho, anchor, n, M, ep); %our
% toc
index = 1;
err_vector = zeros(size(n_vector));
for n = n_vector
    N = n*dim + 1;
    tic
    [U,u, X, XT] = Holger2DEuCallTransform(X_eval,smax, K, T, r, sig1,sig2, rho, anchor, n, M, ep); %our
    toc
    N_spec = 41;
    % Truth
    tic
    True = BSeuCall2D_RBFPUM(X_eval,K,T,r,sig1,sig2,rho,N_spec,M,ep,1,0.15); %Elisabeth
    toc

    trudeau = reshape(True, size(xx));
    UU = reshape(U, size(xx));
    error = (UU- trudeau);
    rel_error = error./trudeau;
    norm_rel_error = max(abs(error(:)));
    err_vector(index) = norm_rel_error
    index = index + 1;
end
q = rate_of_convergence(err_vector, n_vector)
% PLOTTTT!
if exist("XT") == 0
   XT = X; 
end
figure
plot3(X(:,1), X(:,2), u, "ro")
title("Solution at 'Center' Points seen from rotated system")

figure
plot3(XT(:,1), XT(:,2), u, "ro")
title("Solution at 'Center' Points seen from standard system")




figure
hold on
xlim([0.9*Eval_smin, 1.2*Eval_smax])
ylim([0.9*Eval_smin, 1.2*Eval_smax])
plot3(smax*XT(1:n,1),smax*XT(1:n,2),zeros(length(XT(1:n,1))),"r-")
plot3(smax*XT(n+2:N,1),smax*XT(n+2:N,2),zeros(length(XT(1:n,1))),"r-")
f1 = mesh(xx, yy, UU);
title("Solution at Evaluation Points")
view(3)

%True
figure
hold on

mesh(xx, yy, trudeau)
xlim([0.9*Eval_smin, 1.2*Eval_smax])
ylim([0.9*Eval_smin, 1.2*Eval_smax])
% plot3(smax*X(1:n,1),smax*X(1:n,2),zeros(length(X(1:n,1))),"r-")
% plot3(smax*X(n+2:N,1),smax*X(n+2:N,2),zeros(length(X(1:n,1))),"r-")
view(3)
title("True Solution?")

figure
hold on

mesh(xx, yy,error)
xlim([0.9*Eval_smin, 1.2*Eval_smax])
ylim([0.9*Eval_smin, 1.2*Eval_smax])
plot3(smax*XT(1:n,1),smax*XT(1:n,2),zeros(length(XT(1:n,1))),"r-")
plot3(smax*XT(n+2:N,1),smax*XT(n+2:N,2),zeros(length(XT(1:n,1))),"r-")
view(3)
title("Absolute Error")

% figure
% hold on
% mesh(xx, yy,rel_error)
% xlim([0.9*Eval_smin, 1.2*Eval_smax])
% ylim([0.9*Eval_smin, 1.2*Eval_smax])
% plot3(smax*X(1:n,1),smax*X(1:n,2),zeros(length(X(1:n,1))),"r-")
% plot3(smax*X(n+2:N,1),smax*X(n+2:N,2),zeros(length(X(1:n,1))),"r-")
% view(3)
% title("Relative Error")
%% 3D
clear all;
close all;
clc;
%Ecomonic Parameters
K=  20; % Strike
T = 3; % Length of contract
r = 0.02; % Interest rate
sig1 = 0.15;
sig2 = 0.2;
sig3 = 0.1;
rho12 = 0.5; rho13 = 0.5; rho23 = 0.5; 
rho21 = 0.5; rho31 = 0.5; rho32 = 0.5;

C =     [sig1^2,          rho12*sig1*sig2,  rho13*sig1*sig2;
         rho21*sig2*sig1, sig2^2,           rho23*sig1*sig2;
         rho31*sig3*sig1, rho32*sig3*sig2,  sig3^2];

%Numerical Parameters
dim = 2;
n = 20; %Points in each dimention
N = dim*n + 1; %Total nr of points
M = 10; %Number of timesteps.
ep = 10; % Shape parameter
anchor = 0.12; % Anchor, freezing point

% Points for evaluation 
smax = 4*dim*K;        %Largets value for simulation (center points)
Eval_smin = dim*1/3*K;  %Evalutaion min
Eval_smax = dim*5/3*K;  %Evaluation max
temp_x = linspace(Eval_smin,Eval_smax,21);
[xx, yy, zz] = meshgrid(temp_x);
X_eval = [xx(:) yy(:) zz(:)]; %Evaluation points

tic
[U,u, X] = Holger3DEuCall(X_eval, smax, K, T, r, C, anchor, n, M, ep); %our
toc






