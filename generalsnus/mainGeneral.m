clear all;close all; clc;
disp("Saker testas alltid i denna kod, se till att du har kolla på vilka parametrar som du kör med")
%% Ecomonic Parameters
K=  20; % Strike
T = 0.5; % Length of contract
r = 0.02; % Interest rate
sig1 = 0.15; sig2 = 0.2; sig3 = 0.1;
rho12 = 0.5; rho23 = 0.5; rho13 = 0.5;
rho = 0.5;
% 
% Build C-matrix
C = [sig1^2, sig1*sig2*rho12, sig1*sig3*rho13;  %3D
     sig2*sig1*rho12, sig2^2, sig2*sig3*rho23;
     sig3*sig1*rho13, sig3*sig2*rho23, sig3^2];
% % 
% C = [sig1^2, sig1*sig2*rho12; ...         #2D
%     sig2*sig1*rho12, sig2^2];


%% Numerical Parameters

dim = 3;
% maxOrder = min(3,dim-1);
maxOrder = 3;
% C = 0.5*randi(5)*rand(dim);


n = 10; %Points in each dimention
N = dim*n + 1; %Total nr of points Not true :(
M = 10; %Number of timesteps.
ep = 30; % Shape parameter
anchor = 20*ones(dim,1);%[23; 23;23]; % Anchor, freezing point

%% Points for evaluation
smax = 4*dim*K;        %Largets value for simulation (center points)
% Eval_smin = 1/3*K;  %Evalutaion min
% Eval_smax = 5/3*K;  %Evaluation max
% temp_x = linspace(Eval_smin,Eval_smax,11);
% [xx, yy, zz] = meshgrid(temp_x);
% X_eval = [xx(:) yy(:) zz(:)]; %Evaluation points
% X_eval = K*ones(1,dim);

kk = 0;
for k = -20:100
    kk = kk + 1;
    X_eval(kk,:) = (K+k)*ones(1,dim);
end

%% Run
tic
[U, u, X, XT] = generalEuCall(X_eval, smax, K, T, r, C, anchor, n, M, ep, dim, maxOrder);
toc


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
if dim == 22
    N_spec = 41;
    % Truth
    tic
    True = BSeuCall2D_RBFPUM(X_eval,K,T,r,sig1,sig2,rho,N_spec,M,ep,4,0.15); %Elisabeth
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
    plot3(XT(1:n,1),XT(1:n,2),zeros(length(XT(1:n,1))),"r-")
    plot3(XT(n+1:N-1,1),XT(n+1:N-1,2),zeros(length(XT(1:n,1))),"r-")
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
    plot3(XT(1:n,1),XT(1:n,2),zeros(length(XT(1:n,1))),"r-")
    plot3(XT(n+1:N-1,1),XT(n+1:N-1,2),zeros(length(XT(1:n,1))),"r-")
    view(3)
    title("Absolute Error")
end
    %dim plotting
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


