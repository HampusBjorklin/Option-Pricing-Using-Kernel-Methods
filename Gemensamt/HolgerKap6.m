% clear all; close all; clc;
clear all; clc;
n_list = [30];
% n_list = [29];
% n_list = [10, 15, 20, 25, 30];
% n_list = 20:20:100;
n_size = length(n_list);
e = zeros(n_size, 1);
domain = [0, 1];
tic
for n = n_list
    N = 2*n + 1;
    
    anchor = 0.4;
    
    X = zeros(1, N);
    Xgf = zeros(N,1);

    if anchor == domain(1)
        x = linspace(domain(1), domain(2),n+1);
        y = linspace(domain(1), domain(2),n+1);
        x = x(2:end);
        y = y(2:end);
    elseif anchor == domain(2)
        x = linspace(domain(1), domain(2),n+1);
        y = linspace(domain(1), domain(2),n+1);
        x = x(1:end-1);
        y = y(1:end-1);
    else
        pts = round((n+1)*(anchor-domain(1))/(domain(2)-domain(1)));
        l1 = linspace(domain(1), anchor, pts);
        l2 = linspace(anchor, domain(2), n - pts+2);
        x = [l1(1:end-1) l2(2:end)];
        y = [l1(1:end-1) l2(2:end)];
   end

    
    X = [[x',anchor*ones(n,1)]; [anchor, anchor]; [anchor*ones(n,1), y']; [0, 1]; [1,0]];
    
    f = @(x1,x2) x1.^2 + x2.^2 + 13*x2.^4;
%     f = @(x1, x2) basketsolver(x1, x2);
    %f = @(x1,x2) x1 + x2 + cos(x1 + x2);
    % f = @(x1,x2) 10*x1.^2 + 10*x2.^2;
    m = 4;
    for i = 1:N
        Xgf(i) = f(X(i,1), X(i,2));
    end
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
    stepping = (domain(2)-domain(1))/20;  % Om n är en multipel av nämnaren ger det låga fel
    points = [domain(1):stepping:domain(2) ; domain(1):stepping:domain(2)]';
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
%%
convPlot(n_list,e, randi(2000))
figure(randi(2000))
surf(points(:,1), points(:,2), sf-true_val)
title("Error over grid: m = "+ num2str(m) + ", n = " + num2str(n) + ", Anchor = " + num2str(anchor))
figure(randi(2000))
surf(points(:,1), points(:,2), sf)
title("Computed solution")

figure(randi(2000))
surf(points(:,1), points(:,2), true_val)
title("True Solution")


