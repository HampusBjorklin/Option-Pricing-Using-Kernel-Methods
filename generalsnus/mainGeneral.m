clear all;close all; clc;
%% Ecomonic Parameters
K=  20; % Strike
T = 1; % Length of contract
r = 0.02; % Interest rate
sig1 = 0.15; sig2 = 0.2; sig3 = 0.1; sig4 = 0.2; sig5 = 0.12; sig6 = 0.08;
rho12 = 0.5; rho23 = 0.5; rho13 = 0.5;
rho14 = 0.5; rho24 = 0.5; rho34 = 0.5;
rho15 = 0.5; rho25 = 0.5; rho35 = 0.5; rho45 = 0.5;
rho16 = 0.5; rho26 = 0.5; rho36 = 0.5; rho46 = 0.5; rho56 = 0.5;
rho = 0.5;
% 
% Build C-matrix

% 
C = [sig1^2, sig1*sig2*rho12; ...         #2D
    sig2*sig1*rho12, sig2^2];

% C = [sig1^2, sig1*sig2*rho12, sig1*sig3*rho13;  %3D
%      sig2*sig1*rho12, sig2^2, sig2*sig3*rho23;
%      sig3*sig1*rho13, sig3*sig2*rho23, sig3^2];
% % 
% C = [sig1^2, sig1*sig2*rho12, sig1*sig3*rho13, sig1*sig4*rho14;  %4D
%      sig2*sig1*rho12, sig2^2, sig2*sig3*rho23, sig2*sig4*rho24;
%      sig3*sig1*rho13, sig3*sig2*rho23, sig3^2, sig3*sig4*rho34;
%      sig4*sig1*rho13, sig4*sig2*rho24, sig4*sig3*rho34, sig4^2];

% C = [sig1^2, sig1*sig2*rho12, sig1*sig3*rho13, sig1*sig4*rho14,sig1*sig5*rho15 ;  %5D
%      sig2*sig1*rho12, sig2^2, sig2*sig3*rho23, sig2*sig4*rho24,sig2*sig5*rho25;
%      sig3*sig1*rho13, sig3*sig2*rho23, sig3^2, sig3*sig4*rho34, sig3*sig5*rho35;
%      sig4*sig1*rho13, sig4*sig2*rho24, sig4*sig3*rho34, sig4^2, sig4*sig5*rho45;
%      sig5*sig1*rho15, sig5*sig2*rho25, sig5*sig3*rho35, sig5*sig4*rho45, sig5^2];

%  C = [sig1^2, sig1*sig2*rho12, sig1*sig3*rho13, sig1*sig4*rho14,sig1*sig5*rho15, sig1*sig6*rho16;  %6D
%      sig2*sig1*rho12, sig2^2, sig2*sig3*rho23, sig2*sig4*rho24,sig2*sig5*rho25, sig2*sig6*rho26;
%      sig3*sig1*rho13, sig3*sig2*rho23, sig3^2, sig3*sig4*rho34, sig3*sig5*rho35, sig3*sig6*rho36;
%      sig4*sig1*rho13, sig4*sig2*rho24, sig4*sig3*rho34, sig4^2, sig4*sig5*rho45, sig4*sig6*rho46;
%      sig5*sig1*rho15, sig5*sig2*rho25, sig5*sig3*rho35, sig5*sig4*rho45, sig5^2, sig5*sig6*rho56;
%      sig6*sig1*rho16, sig6*sig2*rho26, sig6*sig3*rho36, sig6*sig4*rho46,sig6*sig5*rho56 sig6^2];

%  C = [sig1^2, sig1*sig2*rho12, sig1*sig3*rho13, sig1*sig4*rho14,sig1*sig5*rho15, sig1*sig6*rho16;  %7D
%      sig2*sig1*rho12, sig2^2, sig2*sig3*rho23, sig2*sig4*rho24,sig2*sig5*rho25, sig2*sig6*rho26;
%      sig3*sig1*rho13, sig3*sig2*rho23, sig3^2, sig3*sig4*rho34, sig3*sig5*rho35, sig3*sig6*rho36;
%      sig4*sig1*rho13, sig4*sig2*rho24, sig4*sig3*rho34, sig4^2, sig4*sig5*rho45, sig4*sig6*rho46;
%      sig5*sig1*rho15, sig5*sig2*rho25, sig5*sig3*rho35, sig5*sig4*rho45, sig5^2, sig5*sig6*rho56;
%      sig6*sig1*rho16, sig6*sig2*rho26, sig6*sig3*rho36, sig6*sig4*rho46,sig6*sig5*rho56 sig6^2];


%% Numerical Parameters

dim = 2;
% maxOrder = min(3,dim-1);
maxOrder = 1;
% C = 0.5*ones(dim);


n = 30; %Points in each dimention
N = calcN(dim,maxOrder,n);
disp("Total number of points, N = " + num2str(N));
M = 30; %Number of timesteps.
ep = 50; % Shape parameter
anchor = 20*ones(dim,1);%Anchor, freezing point

%% Points for evaluation
smax = 4*dim*K;        %Largets value for simulation (center points)
Eval_smin = 1/3*K;  %Evalutaion min
Eval_smax = 5/3*K;  %Evaluation max
temp_x = linspace(Eval_smin,Eval_smax,11);
[xx, yy] = ndgrid(temp_x);
X_eval = [xx(:) yy(:)]; %Evaluation points
% X_eval = K*ones(1,dim); % One point. (Removes compute time for evalutation)


%% Run

[U, u, X, XT, time] = generalEuCallOpt(X_eval, smax, K, T, r, C, anchor, n, M, ep, dim, maxOrder);

%% Error and plotting
if dim == 3
    disp("Hallå där, kom ihåg att 'sann' data är gjord med specifika parametrar och är kanske inte samma som dom du kör nu!!")
%     %Find max err.
%     A = load("3DTrue.mat");
%     A = A.Utrue;
%     error = abs(U' - A);
%     disp(max((error)))
    
    figure
    hold on
    plot3(X(:,1),X(:,2),X(:,3),"ro")
    title("Center points (rotated system)");
    view(3)
    
    
    figure
    hold on
    plot3(XT(:,1),XT(:,2),XT(:,3),"ro")
    title("Center points (standard system)");
    plot3(X_eval(:,1),X_eval(:,2),X_eval(:,3),"go")
    view(3)

    figure
    plot(X(1:n,1),u(1:n))
    title("Solution along diagonal direction")

end
if dim == 2
    % Truth
    tic
    True = BSeuCall2D_RBFPUM(X_eval,K,T,r,sig1,sig2,rho,n,M,ep,1,0.15); %Elisabeth
    toc
    
    
    trudeau = reshape(True, size(xx));
    UU = reshape(U, size(xx));
    error = (UU- trudeau);
    err_rel = error ./ norm(trudeau(:),2);
    disp("Maximum absolute error = " + num2str(max(abs(error(:)))));
    disp("Maximum relative error**** = " + num2str(max(abs(err_rel(:)))))
    

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
%     plot3(XT(1:n,1),XT(1:n,2),zeros(length(XT(1:n,1))),"r-")
%     plot3(XT(n+1:N-1,1),XT(n+1:N-1,2),zeros(length(XT(1:n,1))),"r-")
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
%     plot3(XT(1:n,1),XT(1:n,2),zeros(length(XT(1:n,1))),"r-")
%     plot3(XT(n+1:N-1,1),XT(n+1:N-1,2),zeros(length(XT(1:n,1))),"r-")
    view(3)
    title("Absolute Error")
end
    %dim plotting för 1d
    start = 1;
    stop = n;
    x = 1:length(u(1:n));
    for i = 1:dim
        figure
        plot(x,u(start:stop),"-or")
        start = start + n;
        stop = stop + n;
        title("Result in " + num2str(i) + ":th direction")
    end


