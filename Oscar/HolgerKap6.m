% clear all; close all; clc;
clear all; clc;
% n_list = [10:1:50];
n_list = [10];

n_size = length(n_list);
e = zeros(n_size, 1);
domain = [0.00, 1];
Strike = 0.5;

anchor = [0.5, 0.5];
dim = length(domain);
steppingPoints = 11;
tic
for n = n_list
    N = 2*n + 1;
    X = zeros(N, 2);

    startpoints = [[-0.2; anchor(2)], [anchor(1); 1.2]];
    endpoints = [[1.2; anchor(2)], [anchor(1); -0.2]];
%     startpoints = [[-0.2; 0.51], [0.51; 1.2]];
%     endpoints = [[1.2; 0.51], [0.51; -0.2]];
    sz = size(startpoints);
    for i = 1:sz(2)
        for d = 1:dim
            if startpoints(d, i) == endpoints(d, i)
                xi = ones(1, n)*startpoints(d, i);
%             elseif anchor(d) == domain(1)
%                 xi = linspace(domain(1), endpoints(d, i),n+1);
%                 xi = xi(2:end);
%             elseif anchor(d) == domain(2)
%                 xi = linspace(startpoints(d, i), domain(2),n+1);
%                 xi = xi(1:end-1);
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

%     f = @(x1,x2) 10*x1.^2 + 10*x2.^2;
%     Xgf = f(X(:,1), X(:,2));
  
    f = @(x1, x2) basketSolverSingle((x1+x2)/2, (1+x1-x2)/2, Strike);
%     f = @(x1, x2) basketSolverSingle(x1, x2, Strike);
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
    
    
    stepping = (domain(2)-domain(1))/steppingPoints;  % Om n är en multipel av nämnaren ger det låga fel
    points = [domain(1):stepping:domain(2) ; domain(1):stepping:domain(2)]';
    v = @(x,y) [1, -1; 1, 1]/2*([x; y]-[0;1]); 

    [xx, yy] = meshgrid(points(:,1), points(:,2));
    points_x = zeros(size(xx)); points_y = zeros(size(yy));
    for i=1:length(points)
        for j=1:length(points)
            t_xy = v(xx(i, j), yy(i, j));
            points_x(i, j) = t_xy(1);
            points_y(i, j) = t_xy(2);
        end
    end

    sz = length(points(:,1));
    sf = zeros(sz,sz);
    true_val = zeros(sz,sz);
    
    for k1 = 1:length(points(:,1))
        for k2 = 1:length(points(:,2))
            for i = 1:N
                sf(k1,k2)= sf(k1,k2) + alpha(i) * (RepKernel([points_x(k1, k2), points_y(k1, k2)], X(i,:), m, anchor));
            end
            true_val(k1,k2) = f(points_x(k1, k2), points_y(k1, k2));
        end
    end
    
    toc
end
%%

figure(5)
surf(points(:,1), points(:,2), sf-true_val)
% surf((points(:,1)+points(:,2))/2, (1-points(:,1) + flipud(points(:,2)))/2, sf-true_val)
title("Error over grid: m = "+ num2str(m) + ", n = " + num2str(n) + ", Anchor = " + num2str(anchor))
%hold off

figure(4)
%plot3(X(:,1),X(:,2), alpha/max(alpha)*max(sf, [], 'all'),'.')
%hold on
% s2 = surf(0.25 + points(:,1)*0.5, 0.25+points(:,2)*0.5, sf);
s2 = surf(points(:,1), points(:,2), sf);
%hold off
title("Computed solution")

% figure(randi(2000))
figure(7)
surf(points(:,1), points(:,2), true_val)
title("True Solution")

%%
figure(1)
f_xtr = @(x1, x2) [(x1+x2)/2, (1-x1+x2)/2];
xtr = f_xtr(X(:,1), X(:,2));
% xtr = X;
plot(xtr(:,1), xtr(:, 2), '.')
figure(2)
plot3(X(:,1), X(:, 2), Xgf, '.')


[xx, yy] = meshgrid(points(:,1), points(:,2));
figure(3)
for i=1:length(points)
    for j=1:length(points)
        t_xy = f_xtr(xx(i, j), yy(i, j));
        p_x(i, j) = t_xy(1);
        p_y(i, j) = t_xy(2);
    end
end
plot(p_x, p_y, '.')