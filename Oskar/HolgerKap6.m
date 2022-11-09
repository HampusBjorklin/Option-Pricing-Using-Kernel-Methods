clear all; close all; clc;
n_list = [10:1:30];
n_size = length(n_list);
e = zeros(n_size, 1);
for n = n_list
%n = 10; 
N = 2*n + 1;

X = zeros(1, N);
x = linspace(0.05,1,n);
y = linspace(0.05,1,n);
X = [[x',zeros(n,1)]; [0,0]; [zeros(n,1), y']];

f = @(x1,x2) x1.^2 + x2.^2 + 13*x2.^4;
m = 4;
Xgf = f(X(:,1), X(:,2));
A = zeros(N,N);

for i = 1:N
   for j = 1:N
      A(i,j) = RepKernal(X(i,:), X(j,:), m);
   end
end

alpha = A\Xgf;



points = [0:0.1:1 ; 0:0.1:1]';
sz = length(points(:,1));
sf = zeros(sz,sz);
true_val = zeros(sz,sz);


for k1 = 1:length(points(:,1))
    for k2 = 1:length(points(:,2))
        for i = 1:N
            sf(k1,k2)= sf(k1,k2) + alpha(i) * (RepKernal([points(k1,1), points(k2,1)], X(i,:), m));
        end
        true_val(k1,k2) = f(points(k1,1), points(k2,2));
    end
end

e(n==n_list) = norm(true_val - sf, inf); 
disp("Largest error for N = 10^" + num2str(log10(N)) + " is 10^" + num2str(log10(e(n==n_list))));
end

convPlot(n_list,e)