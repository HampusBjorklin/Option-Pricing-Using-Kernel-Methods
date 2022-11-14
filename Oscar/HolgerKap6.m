clear all; close all; clc;
n_list = [10:1:30];
% n_list = [10:15];
% n_list = [10, 15, 20, 25, 30];
% n_list = 20:20:100;
n_size = length(n_list);
e = zeros(n_size, 1);
domain = [0, 100];
tic
for n = n_list
    %n = 10; 
    N = 2*n + 1;
    
    
    anchor = 50;
    
    X = zeros(1, N);
    x = linspace(domain(1), domain(2),n+1);
    y = linspace(domain(1), domain(2),n+1);
    
    anchInd = find(x >= anchor, 1);
    % anchor = x(anchInd);
    
    if anchor == domain(1)
        x = x(2:end);
        y = y(2:end);
    elseif anchor == domain(2)
        x = x(1:end-1);
        y = y(1:end-1);
    else
        if min(size(find(x >= anchor, 1))) == 0
            pts = round((n+1)*(anchor-domain(1))/(domain(2)-domain(1)));
            l1 = linspace(domain(1), anchor, pts);
            l2 = linspace(anchor, domain(2), n - pts + 1);
            x = [l1 l2(2:end)];
            y = [l1 l2(2:end)];
%             x = [x(1:anchInd-1), anchor, x(anchInd: end)];
%             y = [y(1:anchInd-1), anchor, y(anchInd: end)];

        else
            x = [x(1:anchInd-1), x(anchInd+1: end)];
            y = [y(1:anchInd-1), y(anchInd+1: end)];
        end
    end
    
    X = [[x',anchor*ones(n,1)]; [anchor, anchor]; [anchor*ones(n,1), y']];
    
    % f = @(x1,x2) x1.^2 + x2.^2 + 13*x2.^4;
    % f = @(x1,x2) x1 + x2 + cos(x1 + x2);
    f = @(x1,x2) 10*x1.^2 + 10*x2.^2;
    m = 4;
    Xgf = f(X(:,1), X(:,2));
    % Xgf = basketsolver(X);
    % f = @(x1, x2) basketSolverSingle([x1, x2]);
    A = zeros(N,N);
    
    for i = 1:N
       for j = 1:N
          A(i,j) = RepKernel(X(i,:), X(j,:), m, anchor);
       end
    end
    
    alpha = A\Xgf;
    
    
    % points = [domain(1):0.1:domain(2) ; domain(1):0.1:domain(2)]';
    points = [domain(1):5:domain(2) ; domain(1):5:domain(2)]';
    sz = length(points(:,1));
    sf = zeros(sz,sz);
    true_val = zeros(sz,sz);
    
    
    
    for k1 = 1:length(points(:,1))
        for k2 = 1:length(points(:,2))
            for i = 1:N
                sf(k1,k2)= sf(k1,k2) + alpha(i) * (RepKernel([points(k1,1), points(k2,1)], X(i,:), m, anchor));
            end
            % sf(k1,k2)= sum(alpha .* (RepKernel([points(k1,1), points(k2,1)], X(:,:), m, anchor)));
            true_val(k1,k2) = f(points(k1,1), points(k2,2));
        end
    end
    
    % e(n==n_list) = norm(true_val - sf, inf);
    e(n==n_list) = max(max(abs(true_val - sf)));
    disp("Largest error for N = 10^" + num2str(log10(N)) + " is 10^" + num2str(log10(e(n==n_list))));
    disp(max(max(abs(true_val - sf))))
    toc
end

convPlot(n_list,e)