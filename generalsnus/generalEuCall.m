function [U, u, X, XT] = generalEuCall(X_eval, smax, K, T, r, C, anchor, n, M, ep, dim, maxOrder)
% Calculates the fair option price at points X_eval using a reproducing
% kernel and radial basis funtions. This one is in 3D :)


% Kolla på Coeff functionerna och RepKernal
%Check parameter matrix
if issymmetric(C) == 0 || length(C) ~= dim
    error("Parameter matrix is not symmetric, please fix :)")
end

Ne = length(X_eval(:,1));    %Number of evaluation points
N = dim*n + 1;          %Total number of centerpoints


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

% tm1 = linspace(0,1,n);
% tm2 = linspace(-0.5,0.5,n);
% [xx1,yy1] = meshgrid(tm1,tm2); % centers of partitions

% X = [xx1(:) yy1(:)];
% plot(X(:,1),X(:,2),"go")

N = length(X);
XT = f_v2s(X);              %In s

X_eval = f_s2v(X_eval);     %In v


% Define boundary points
distClose = 0.01;
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

% Coeff -> Nodal (Always one)
%First Order dervs A_1, A_2,...,A_d
[firstOrder(1:NFirstOrderTerms).matrix] = deal(zeros(N));
%Second order dervs A_11, A_22,.. ,A12,...., Ad-1d
[secondOrder(1:NSecondOrderTerms).matrix] = deal(zeros(N));


%% Generate the matrices
%A0

[xx, yy] = meshgrid(1:N);
A0 = GeneralRepKernel(X(xx(:),:), X(yy(:),:), ep, maxOrder);
A0 = (reshape(A0,N,N))';

%Generate all derivative matrices. Note that there might be more mixed
%derivatives than pure derivatives. Hence the split in for. The
%secoundOrder struct is also structured followingly:
% First "NPureSecondOrderTerms"- Pure derivatives then "NMixedTerms" mixed
% derivative (this explains indexing)
for k = 1:NFirstOrderTerms
    tmp = GeneralRepKernelFirstDer(X(xx(:),:), X(yy(:),:), ep, k, maxOrder);
    firstOrder(k).matrix = (reshape(tmp,N,N))';
    tmp = GeneralRepKernelSecondDer(X(xx(:),:), X(yy(:),:), ep, k, maxOrder);
    secondOrder(k).matrix = (reshape(tmp,N,N))';
end
for k = 1:NMixedTerms
    tmp = GeneralRepKernelSecondMixed(X(xx(:),:), X(yy(:),:), ep, mixers(k,:), maxOrder);
    secondOrder(k + NPureSecondOrderTerms).matrix = (reshape(tmp,N,N))';
end

%% Interiour Operators with Transform system
%These have size R^(NInter x N)
B0 = A0(indInter, :);
%Allocate
[BMatrices(1:(NFirstOrderTerms + NSecondOrderTerms)).matrix] = deal(zeros(N));

% The general derivative is a linear combination of all derivatives in same
% order: ds1 = [c1 c2 c3 c4] [dv1 dv2 dv3 dv4]'
% Note that derivatives are matrices so special "vecmult" is required. They
% are however constant for all s-derivatives. Only thing changing are the
% "c".

%Allocate
D1_mat = zeros(NInter,N*NFirstOrderTerms);
D2_mat = zeros(NInter,N*NSecondOrderTerms);

%The elements in the coeff-vector is taken from transformation matrix
%Eventuellt så är detta lite oklart, men jävligt snygg one-liners
% R^(1 x NFirstOrderTerms)
D1_coeff = @(i) [MS2V(1:dim,i)]';
% R^(1 x NSecondOrderTerms)
D2_coeff = @(i,j) [(MS2V(1:dim,i).*MS2V(1:dim,j)); MS2V(mixers(:,1),i).*MS2V(mixers(:,2),j) + MS2V(mixers(:,2), i).*MS2V(mixers(:,1), j)]';

%Fill the Derivative Vector/Matrix
for i = 1:NFirstOrderTerms
    D1_mat(:,1 + (i-1)*N:i*N) = [firstOrder(i).matrix(indInter,:)];
end
for i = 1:NSecondOrderTerms
    D2_mat(:,1 + (i-1)*N:i*N) =  [secondOrder(i).matrix(indInter,:)];
end


%The B Matrices (S-derivatives in terms of V) are built
for i = 1:NFirstOrderTerms
    BMatrices(i).matrix = special_matmult(D1_coeff(i), D1_mat);
    BMatrices(i + NFirstOrderTerms).matrix = special_matmult(D2_coeff(i,i), D2_mat);
end
for i = 1:NMixedTerms
    BMatrices(i + NFirstOrderTerms + NPureSecondOrderTerms).matrix ...
        = special_matmult(D2_coeff(mixers(i,1),mixers(i,2)), D2_mat);
end

%Diagonal S-Matrices.
[SMatrices(1:dim).matrix] = deal(zeros(NInter,NInter));
for i = 1:dim
    SMatrices(i).matrix = diag(XT(indInter,i));
end

%% Build the BS-Operator
% dim terms are +r*V_i*D_i                      (1)
% dim terms are +0.5*V_i^2*C(i,i)*D_ii          (2)
% NMixed terms are    +V_i*V_j*C(i,j)*D_ij      (3)
% 1 term -r*D0                                  (4)
Operator = zeros(NInter,N);
for i = 1:dim               %(1) (2)
    Operator = (Operator + r*SMatrices(i).matrix*BMatrices(i).matrix);
    Operator = Operator + 0.5*C(i,i)*(SMatrices(i).matrix).^2 ...
        * BMatrices(i + NFirstOrderTerms).matrix;
end
for i = 1:NMixedTerms %(3)
    ii = mixers(i,1);
    jj = mixers(i,2);
    Operator = Operator + C(ii,jj)*SMatrices(ii).matrix*SMatrices(jj).matrix ...
        *BMatrices(i + NFirstOrderTerms + NPureSecondOrderTerms).matrix;
end
Operator = (Operator - r*B0)/A0; %(4)

%Add Interiour operator at correct indices in the Large operator
L = zeros(N,N);
L(indInter, :) = L(indInter, :) + Operator;


%% Time Solver
% Using the BDF2 function.

% k = T/(M-1);
% beta0 = k* 2/3;
% beta1 = 4/3;
% beta2 = 1/3;
[k,beta0,beta1,beta2]=BDF2coeffs(T,M);
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
    rhs = beta1(nextstep)*u - beta2(nextstep)*u0;
    
    tn = tvec(nextstep); % Should be one step ahead
    
    rhs(indClose,:) = closeBC(XT(indClose,:),K,r,tn);
    rhs(indFar,:) = farBC(XT(indFar,:),K,r,tn);
    
    % move the solutions one step
    u0 = u;
end
end

%% Obtain solution at Desired Evalutation Points
% Evaluation matrix

[xx, yy] = meshgrid(1:Ne,1:N);
E = GeneralRepKernel(X_eval(xx(:),:), X(yy(:),:), ep, maxOrder);
E = reshape(E,N,Ne)'; %Note that this wierdness is due to reshape placing columns first (we want rows :))
                      % Not problem since other matrices are symemetric. lo
E = E/A0;
% Solution
U = E*u;
%Rescale
X_eval = X_eval*smax;
U = U*smax; K = K*smax;
X = X.*smax; XT = XT.*smax;
u = u*smax;
end