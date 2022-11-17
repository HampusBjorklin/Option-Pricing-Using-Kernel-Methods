% clear all; close all; clc;
clear all; clc;
% n_list = [10:1:50];
n_list = [30];
% n_list = [10, 15, 20, 25, 30];
% n_list = 20:20:100;
n_size = length(n_list);
e = zeros(n_size, 1);
domain = [0.00, 1];
Strike = 0.5;
% anchor = 0.5;
% anchor = [0.5, 0.5];
anchor = [0.55, 0.55];
% anchor = [0.9, 0.9];
dim = length(domain);
steppingPoints = 11;
tic
for n = n_list
    %n = 10; 
    % N = 2*n + 1;
    N = 3*n + 1;
    % N = 4*n + 1;
    % N = N + 2;
    X = zeros(N, 2);
%     startpoints = [[0, anchor(2)]; [anchor(1), 0]];  % storlek linjer*dimension
%     endpoints = [[1, anchor(2)]; [anchor(1), 1]];
    theta = deg2rad(45);
    startpoints = rotatePoint([[anchor(1), 0]; [1, anchor(2)]]', anchor', theta)';  % storlek linjer*dimension
    endpoints = rotatePoint([[anchor(1), 1]; [0, anchor(2)]]', anchor', theta)';
%     startpoints = rotatePoint([[0.5, 0]; [1, 0.5]]', anchor', theta)';  % storlek linjer*dimension
%     endpoints = rotatePoint([[0.5, 1]; [0, 0.5]]', anchor', theta)';
    startpoints = [startpoints, [anchor(1);0]];
    endpoints = [endpoints, [anchor(1);1]];
%     anchor = [anchor anchor(2)];
%     startpoints = [[0, 0]; [1, 0]];  % storlek linjer*dimension
%     endpoints = [[1, 1]; [0, 1]];
%     startpoints = [[0, 0]; [1, 0.8]];  % storlek linjer*dimension
%     endpoints = [[1, 1]; [0.8, 1]];
    sz = size(startpoints);
    for i = 1:sz(2) % byt till rätt sätt att iterera
        for d = 1:dim
            if startpoints(d, i) == endpoints(d, i)
                xi = ones(1, n)*startpoints(d, i);
            elseif anchor(d) == domain(1)
                xi = linspace(domain(1), endpoints(d, i),n+1);
                xi = xi(2:end);
            elseif anchor(d) == domain(2)
                xi = linspace(startpoints(d, i), domain(2),n+1);
                xi = xi(1:end-1);
            else
                pts = round(abs((n+1)*(anchor(d)-startpoints(d, i)))...
                    /abs((endpoints(d, i)-startpoints(d, i))));
                l1 = linspace(startpoints(d, i), anchor(d), pts);
                l2 = linspace(anchor(d), endpoints(d, i), n - pts+2);
                xi = [l1(1:end-1) l2(2:end)];
            end
            X(1+(i-1)*n:i*n, d) = xi;
        end
    end
    X(end, :) = anchor;

%     f = @(x1,x2) x1.^2 + x2.^2 + 13*x2.^4;
%     f = @(x1,x2) x1.^2 + x2.^2 + 13*x2.^4 + 0*((2)).*(x1>x2);
%     f = @(x1,x2) x1 + x2 + cos(x1 + x2);
%     f = @(x1,x2) x1.^2-x2.^2 + (x1>0.5).*ones(size(x1)); 
% 
%     f = @(x1,x2) 10*x1.^2 + 10*x2.^2;
%     Xgf = f(X(:,1), X(:,2));
    
    f = @(x1, x2) basketSolverSingle(x1, x2, Strike);
    Xgf = zeros(N, 1);
    for i=1:N
        Xgf(i) = f(X(i, 1), X(i, 2));
    end

    m = 1;
    A = zeros(N,N);
    
    for i = 1:N
       for j = 1:N
          A(i,j) = RepKernel(X(i,:), X(j,:), m, anchor);
       end
    end
    
    alpha = A\Xgf;
    
    
    % points = [domain(1):0.1:domain(2) ; domain(1):0.1:domain(2)]';
    stepping = (domain(2)-domain(1))/steppingPoints;  % Om n är en multipel av nämnaren ger det låga fel
    points = [domain(1):stepping:domain(2) ; domain(1):stepping:domain(2)]';
    sz = length(points(:,1));
    sf = zeros(sz,sz);
    true_val = zeros(sz,sz);
    
    
    
    for k1 = 1:length(points(:,1))
        for k2 = 1:length(points(:,2))
            for i = 1:N
                sf(k1,k2)= sf(k1,k2) + alpha(i) * (RepKernel([points(k1,1), points(k2,2)], X(i,:), m, anchor));
            end
            % sf(k1,k2)= sum(alpha .* (RepKernel([points(k1,1), points(k2,1)], X(:,:), m, anchor)));
            true_val(k1,k2) = f(points(k1,1), points(k2,2));
        end
    end
    
    % e(n==n_list) = norm(true_val - sf, inf);
%     e(n==n_list) = max(max(abs(true_val - sf)));
%     disp("Largest error for N = 10^" + num2str(log10(N)) + " is 10^" + num2str(log10(e(n==n_list))));
%     disp(max(max(abs(true_val - sf))))
    toc
end
%%
% % convPlot(n_list,e, randi(2000))
% figure(randi(2000))
% %plot3(repmat([1, 1]/sqrt(2), [N,1])*X, repmat([1, -1]/sqrt(2), [N,1])*X, alpha/max(alpha)*max(sf-true_val, [], 'all'),'.')
% %hold on
figure(5)
surf(points(:,1), points(:,2), sf-true_val)
title("Error over grid: m = "+ num2str(m) + ", n = " + num2str(n) + ", Anchor = " + num2str(anchor))
%hold off

figure(4)
%plot3(X(:,1),X(:,2), alpha/max(alpha)*max(sf, [], 'all'),'.')
%hold on
surf(points(:,1), points(:,2), sf)
%hold off
title("Computed solution")

% figure(randi(2000))
figure(7)
surf(points(:,1), points(:,2), true_val)
title("True Solution")

%%
figure(1)
plot(X(:,1), X(:, 2), '.')
figure(2)
plot3(X(:,1), X(:, 2), Xgf, '.')