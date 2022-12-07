
Nx = 20;
Ny = 20;
xxx = linspace(0,1,100)
xx = linspace(0,1,Nx);
yy = linspace(0,1,Ny);

[X, Y] = meshgrid(xx, yy)



mq_rbf = @(x,y,eps) sqrt(1 + (eps*norm(x-y,2)));

% test_function = @(x) sin(2*pi*x);
test_function = @(x) max((x-0.5), 0);

A = zeros(length(xx), length(yy));
eps = 0.5;

for i = 1:length(xx)
    for j = 1:length(yy)
        A(i, j) = mq_rbf(xx(i), yy(j), eps);
    end
end

disp(A)

u = test_function(xx)';
lambda = A\u;

J = zeros(1, Nx);
for i = 1:Nx
    sum = 0;
    for j = 1:Nx
        sum = sum + lambda(j)*mq_rbf(xx(i), xx(j), eps);
    end
    J(i) = sum
end

plot(xx, J)
hold on
plot(xxx, test_function(xxx))