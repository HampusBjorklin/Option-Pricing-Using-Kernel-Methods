function [U, u, X, XT] = generalEuCall(X_eval, smax, K, T, r, C, anchor, n, M, ep, dim, maxOrder)
% Calculates the fair option price at points X_eval using a reproducing
% kernel and radial basis funtions
tot_time = 0;
%Check parameter matrix
if issymmetric(C) == 0 || length(C) ~= dim
    error("Parameter matrix is not symmetric, please fix :)")
end

Ne = length(X_eval(:,1));    %Number of evaluation points
N = calcN(dim,maxOrder,n);          %Total number of centerpoints


% Scale down the Evalutaion points to region [0, 1]
X_eval = X_eval./smax;
K = K /smax;
anchor = anchor/smax;



% %Transform matrices & functions
MV2S = GCV2S(dim);
MS2V = MV2S';


f_s2v = @(S) (MS2V*S')';
f_v2s = @(V) (MV2S*V')';

%Find in V
anchorT = MS2V*anchor;


%Center points
X = getXVector(anchorT,n,maxOrder);  %In v

XT = f_v2s(X);              %In s

X_eval = f_s2v(X_eval);     %In v


% Define boundary points
distClose = 0.01; %Computers cant count
distFar = 0.99;
range = sqrt(sum(XT.^2,2));             
indClose = find(range<=distClose);
indFar = find(range>=distFar);
indInter = find(range>distClose & range<distFar);


% Boundary Conditions
closeBC = @(S,K,r,t) zeros(size(S,1),1);
farBC =  @(S,K,r,t) max(1/dim * sum(S,2) - K*exp(-r*t),0); %Call on me :)


%Find interior points
XInter = X(indInter,:);
NInter = length(indInter);


%% Differentiation matrices

%General number of terms
NFirstOrderTerms = dim;
NPureSecondOrderTerms = dim;
NMixedTerms = nchoosek(dim,2);
NSecondOrderTerms = dim + NMixedTerms;

mixers = nchoosek(1:dim,2); %Create the order of mixed terms
all_d2 = [[1:dim]',[1:dim]';mixers];

%% Rotation matris
%The elements in the coeff-vector is taken from transformation matrix
%Eventuellt så är detta lite oklart, men jävligt snygg one-liners
% R^(1 x NFirstOrderTerms)
D1_coeff = @(derDim) [MS2V(1:dim,derDim)]';
% R^(1 x NSecondOrderTerms)
D2_coeff = @(derDim) [(MS2V(1:dim,derDim(:,1)).*MS2V(1:dim,derDim(:,2))); MS2V(mixers(:,1),derDim(:,1)).*MS2V(mixers(:,2),derDim(:,2)) + MS2V(mixers(:,2), derDim(:,1)).*MS2V(mixers(:,1), derDim(:,2))]';

% [B1; B2; ...; B_d] = D1_coeff_mat * [A1; A2; ...;]
D1_coeff_mat = D1_coeff(1:dim);
D2_coeff_mat = D2_coeff(all_d2);

%% Här testas det!
tic
[xx, yy] = meshgrid(1:N);
A0 = GeneralRepKernel(X(xx(:),:), X(yy(:),:), ep, maxOrder);
A0 = (reshape(A0,N,N))';
[xx, yy] = meshgrid(indInter,1:N);
time = toc;
disp("Generating A0 took " + num2str(time) + "s")
tot_time = tot_time + time;
tic
Operator = - r*A0(indInter, :);
%Diagonal S-Matrices.
S1 = XT(indInter,:);
S2 = [XT(indInter,:).^2,XT(indInter,mixers(:,1)).*XT(indInter,mixers(:,2))];
sig_D2 = [0.5*diag(C);zeros(NMixedTerms,1)];
for i = 1:NMixedTerms
   sig_D2(i + NPureSecondOrderTerms) = C(mixers(i,1), mixers(i,2));
end

%Allocate
tic
for i = 1:NFirstOrderTerms
    tmp = (reshape(GeneralRepKernelFirstDer(X(xx(:),:), X(yy(:),:), ep, i, maxOrder),N,NInter))';
    Operator = Operator + r*(S1*D1_coeff_mat(:,i)).*tmp;
end
time = toc;
disp("Generating D1 of Operator took " + num2str(time) + "s")
tot_time = tot_time + time;
tic 
for i = 1:NSecondOrderTerms
    if i <= NPureSecondOrderTerms
        tmp = (reshape(GeneralRepKernelSecondDer(X(xx(:),:), X(yy(:),:), ep, i, maxOrder),N,NInter))';
    else
        tmp = (reshape(GeneralRepKernelSecondMixed(X(xx(:),:), X(yy(:),:), ep, mixers(i - NPureSecondOrderTerms,:), maxOrder),N,NInter))';
    end
    Operator = Operator + (S2*(D2_coeff_mat(:,i).*sig_D2)).*tmp;
end
time = toc;
disp("Generating D2 of Operator took " + num2str(time) + "s")
tot_time = tot_time + time;
tic 

Operator = (Operator)/A0;
time = toc;
disp("Turning Operator to nodal rep. took " + num2str(time) + "s")
tot_time = tot_time + time;

%Add Interiour operator at correct indices in the Large operator
L = zeros(N,N);
L(indInter, :) = L(indInter, :) + Operator;


%% Time Solver
% Using the BDF2 function.
tic
k = T/(M-1);
beta0 = k* 2/3;
beta1 = 4/3;
beta2 = 1/3;

%Used to elimitate BC
C = eye(N) -beta0*L;

%Initialization
u0 = farBC(XT, K, r, 0); %IC
u = u0;
rhs = u0;
tvec = 0:k:T;

[Lc,Uc] = lu(C); %Split Matrix to increase speeeeeed

%Time steppin'
if T ~= 0
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
end

time = toc;
disp("Time solver took "+ num2str(time) + "s")
tot_time = tot_time + time;

%% Obtain solution at Desired Evalutation Points
% Evaluation matrix
tic
[xx, yy] = meshgrid(1:Ne,1:N);
E = GeneralRepKernel(X_eval(xx(:),:), X(yy(:),:), ep, maxOrder);
E = reshape(E,N,Ne)'; %Note that this wierdness is due to reshape placing columns first (we want rows :))

E = E/A0;
time = toc;
disp("Evaluation took " + num2str(time) + "s")
% Solution
U = E*u;
%Rescale back
X_eval = X_eval*smax;
U = U*smax; K = K*smax;
X = X.*smax; XT = XT.*smax;
u = u*smax;
disp("Total run time: " + num2str(tot_time) + "s")
end