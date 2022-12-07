eps = 0.5;
Nx = 10;
Ny = 10;
xxx = linspace(-2, 2,100);
yyy = linspace(-2, 2, 100);

xx = linspace(-2,2,Nx);
yy = linspace(-2,2,Ny);

[X1, Y1] = meshgrid(xxx, yyy);
[X2, Y2] = meshgrid(xx, yy);
nodes = zeros(Nx*Ny, 2);
for i = 1:Nx
    for j = 1:Ny
        nodes((i-1)*Ny + j,:) = [xx(i), yy(j)];
    end
end

nodes

mq_rbf = @(x,y,eps) sqrt(1 + (eps*norm(x-y,2)));
test_function = @(x, y) x.*exp(-x.^2-y.^2);

F1 = test_function(X1, Y1);
surf(X1, Y1, F1)
hold on
u = zeros(length(nodes),1);
A = zeros(length(nodes), length(nodes));

for i = 1:length(nodes)
    for j = 1:length(nodes)
        A(i, j) = mq_rbf(nodes(i,:), nodes(j,:), eps);
    end
end

for i = 1:length(nodes)
    u(i) = test_function(nodes(i,1), nodes(i,2));
end

lambda = A\u;
J = zeros(Nx, 1);

for i = 1:length(nodes)
    sum = 0;
    for j = 1:length(nodes)
        sum = sum + lambda(j)*mq_rbf(nodes(i,:), nodes(j,:), eps);
    end
    J(i) = sum;
end
disp(J)
U = zeros(Nx, Ny);
for i = 1:Nx
    for j = 1:Ny
        U(j, i) = J((i-1)*Ny+j);
    end
end
surf(X2, Y2, U)

