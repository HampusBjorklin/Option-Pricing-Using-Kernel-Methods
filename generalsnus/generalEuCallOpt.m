function [U, u, X, XT, tot_time] = generalEuCallOpt(X_eval, smax, K, T, r, C, anchor, n, M, ep, dim, maxOrder)
% Calculates the fair option price at points X_eval using a reproducing
% kernel and radial basis funtions
Ntime2understand = 10;

tot_time = 0;
eps2 = ep^2;
dim_list = 1:dim;
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
if dim == maxOrder
   MV2S = eye(dim); 
end

MS2V = MV2S';


f_s2v = @(S) (MS2V*S')';
f_v2s = @(V) (MV2S*V')';


anchorT = MS2V*anchor;

limit = [zeros(dim,1),ones(dim,1)];
if dim ~= maxOrder
for i = 2:dim
    v = MV2S(:,i);
    lim = min(abs(anchor./v));
    limit(i,:) = [-lim, lim];
end
end
%Center points
X = getXVector(anchorT,n,maxOrder, limit);  %In v

XT = f_v2s(X);              %In s
X_eval = f_s2v(X_eval);     %In v


%Drop the points that are negativ
indices = find(XT < -0.00001); % find the indices of negative elements
[row, col] = ind2sub(size(XT), indices); % convert indices to row-column form
unique_rows = unique(row); % find the unique rows containing negative elements
X(unique_rows,:) = [];
XT(unique_rows,:) = [];
for j = 1:Ntime2understand
disp("Dropped " +num2str(length(unique_rows)) + " negative elements")
end
N = length(XT);


% Define boundary points
distClose = 0.01; %Computers cant count
distFar = 0.99;
range = sqrt(sum(XT.^2,2));
indClose = find(range<=distClose);
indFar = find(range>=distFar);
indInter = find(range>distClose & range<distFar);

% indFar = [find(range>=distFar); (n+1:n:n*dim)'; (2*n:n:n*dim)'];
% [~, indInter] = setdiff(XT, [XT(indFar, :); XT(indClose, :)], 'rows');

% Boundary Conditions
closeBC = @(S,K,r,t) zeros(size(S,1),1);
farBC =  @(S,K,r,t) max(1/dim * sum(S,2) - K*exp(-r*t),0); %Call on me :)


%Find interior points
XInter = X(indInter,:);
NInter = length(indInter);


%% 
%General number of terms
NFirstOrderTerms = dim;
NPureSecondOrderTerms = dim;
NMixedTerms = nchoosek(dim,2);
NSecondOrderTerms = dim + NMixedTerms;

mixers = nchoosek(1:dim,2);   %Create the order of mixed terms
all_d1 = [1:dim];
all_d2 = [[1:dim]',[1:dim]';mixers];%All 2:nd derivatives in order

%% Derivative matrix
%The matrix transforming the derivatives in s to derivative in v
% [B1; B2; ...; B_d] = D1_coeff_mat * [A1; A2; ...;]

%Functions for each element
D1_coeff_func = @(derDim) [MS2V(1:dim,derDim)]';
D2_coeff_func = @(derDim) [(MS2V(1:dim,derDim(:,1)).*MS2V(1:dim,derDim(:,2))); MS2V(mixers(:,1),derDim(:,1)).*MS2V(mixers(:,2),derDim(:,2)) + MS2V(mixers(:,2), derDim(:,1)).*MS2V(mixers(:,1), derDim(:,2))]';
%Matrices
D1_coeff_mat = D1_coeff_func(all_d1);
D2_coeff_mat = D2_coeff_func(all_d2);

%% Building the operator (With less memory)
% With reproducing kernal and multiquadric kernal

%Multiquadric reproducing kernel function
% multi = @(a,b) sqrt(1+eps2*(a-b).^2);


%Elements in the NxN matrix
tic
[xx, yy] = meshgrid(1:N);
x_long = X(xx(:),:);
y_long = X(yy(:),:);

%A0
A0 = ones(size(x_long,1), 1); %+1 is always included

s = {};

for i=1:maxOrder
    subsets = nchoosek(dim_list,i);
    for j = 1:size(subsets,1)
        s(end+1) = {subsets(j,:)};
    end
end
for i=1:length(s)
    arr = cell2mat(s(i));
    arr_len = length(arr);
    coeff = 1;
    for j=1:arr_len
%         coeff = coeff.*multi(x_long(:,arr(j)),y_long(:,arr(j)));
        coeff = coeff.*(1+eps2*(x_long(:,arr(j))-y_long(:,arr(j))).^2);
    end
    A0 = A0 + sqrt(coeff);
end

A0 = (reshape(A0,N,N))';
[xx, yy] = meshgrid(indInter,1:N);
x_long = X(xx(:),:);
y_long = X(yy(:),:);

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

tic
for ii = 1:NFirstOrderTerms
    derivative_coeff = (eps2)*(x_long(:,dim_list(ii))-y_long(:,dim_list(ii)))./(sqrt(eps2*(x_long(:,dim_list(ii))-y_long(:,dim_list(ii))).^2+1));
    s = {};
    d = dim_list;
    d=d(d~=dim_list(ii));
    for i=1:maxOrder-1
        subsets = nchoosek(d,i);
        for j = 1:size(subsets,1)
            s(end+1) = {subsets(j,:)};
        end
    end
    tmp = ones(size(x_long,1), 1); %+1 is always included
    for i=1:length(s)
        arr = cell2mat(s(i));
        arr_len = length(arr);
        coeff = 1;
        for j=1:arr_len
%             coeff = coeff.*multi(x_long(:,arr(j)),y_long(:,arr(j)));
              coeff = coeff.*(1+eps2*(x_long(:,arr(j))-y_long(:,arr(j))).^2);
        end
        tmp = tmp + sqrt(coeff);
    end
    tmp = derivative_coeff.*tmp;
    tmp = (reshape(tmp,N,NInter))';
    Operator = Operator + r*(S1*D1_coeff_mat(:,ii)).*tmp;
end

time = toc;
disp("Generating D1 of Operator took " + num2str(time) + "s")
tot_time = tot_time + time;
tic

%2nd Order terms for operator
for ii = 1:NSecondOrderTerms
    if ii <= NPureSecondOrderTerms
        derivative_coeff = (eps2)./((eps2*(x_long(:,dim_list(ii))-y_long(:,dim_list(ii))).^2+1).^(3/2));
        s = {};
        d = dim_list; %Vector with all dimensions except for the one we look at
        d=d(d~=dim_list(ii));
        
        for i=1:maxOrder-1
            subsets = nchoosek(d,i);
            for j = 1:size(subsets,1)
                s(end+1) = {subsets(j,:)};
            end
        end
        tmp = ones(size(x_long,1), 1); %+1 is always included
        for i=1:length(s)
            arr = cell2mat(s(i));
            arr_len = length(arr);
            coeff = 1;
            for j=1:arr_len
%                 coeff = coeff.*multi(x_long(:,arr(j)),y_long(:,arr(j)));
                coeff = coeff.*(1+eps2*(x_long(:,arr(j))-y_long(:,arr(j))).^2);
            end
            tmp = tmp + sqrt(coeff);
        end
        tmp = derivative_coeff.*tmp;
        tmp = (reshape(tmp,N,NInter))';
    else
        derivative_coeff = (eps2^2).*(x_long(:,mixers(ii - NPureSecondOrderTerms,1))-y_long(:,mixers(ii - NPureSecondOrderTerms,1)))./(sqrt(eps2*(x_long(:,mixers(ii- NPureSecondOrderTerms,1))-y_long(:,mixers(ii- NPureSecondOrderTerms,1))).^2+1)) .*...
            (x_long(:,mixers(ii- NPureSecondOrderTerms,2))-y_long(:,mixers(ii- NPureSecondOrderTerms,2)))./(sqrt(eps2.*(x_long(:,mixers(ii- NPureSecondOrderTerms,2))-y_long(:,mixers(ii- NPureSecondOrderTerms,2))).^2+1));

        s = {};
        d = dim_list; %Vector with all dimensions except for the one we look at
        d=d(d~=mixers(ii- NPureSecondOrderTerms,1));
        d=d(d~=mixers(ii- NPureSecondOrderTerms,2));
        for i=1:maxOrder-2
            subsets = nchoosek(d,i);
            for j = 1:size(subsets,1)
                s(end+1) = {subsets(j,:)};
            end
        end
        if isempty(s) && maxOrder<2
            tmp = zeros(size(x_long,1), 1);
        else
            tmp = ones(size(x_long,1), 1); %+1 is always included
            for i=1:length(s)
                arr = cell2mat(s(i));
                arr_len = length(arr);
                coeff = 1;
                for j=1:arr_len
%                     coeff = coeff.*multi(x_long(:,arr(j)),y_long(:,arr(j)));  
                    coeff = coeff.*(1+eps2*(x_long(:,arr(j))-y_long(:,arr(j))).^2);
                end
                tmp = tmp + sqrt(coeff);
            end
            tmp = derivative_coeff.*tmp;
        end
        tmp = (reshape(tmp,N,NInter))';
    end
    Operator = Operator + (S2*(D2_coeff_mat(:,ii).*sig_D2)).*tmp;
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
% [k, beta0, beta1, beta2] = BDF2coeffs(T,M);
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
x_long = X_eval(xx(:),:);
y_long = X(yy(:),:);

s = {};

d = 1:dim; %Vector with all dimensions
for i=1:maxOrder
    subsets = nchoosek(d,i);
    for j = 1:size(subsets,1)
        s(end+1) = {subsets(j,:)};
    end
end

E = ones(size(x_long,1), 1); %+1 is always included

for i=1:length(s)
    arr = cell2mat(s(i));
    arr_len = length(arr);
    coeff = 1;
    for j=1:arr_len
        coeff = coeff.*(1+eps2*(x_long(:,arr(j))-y_long(:,arr(j))).^2);
    end
    E = E + sqrt(coeff);
end
E = reshape(E,N,Ne)'; %

E = E/A0;
time = toc;
disp("Evaluation took " + num2str(time) + "s")
% Solution
U = E*u;

%Rescale
X_eval = X_eval*smax;
U = U*smax; K = K*smax;
X = X.*smax; XT = XT.*smax;
u = u*smax;
disp("Total run time: " + num2str(tot_time) + "s")
end