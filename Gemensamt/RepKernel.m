function [K] = RepKernel(x, y , m, anchor)
% Compute the repuducing kernal for Lambda = {null, {1}, {2}}
% X = [x1,..,xn]; y = [y1,...y]
    K = 1 + kernel(x(1), y(1), m, anchor) +  kernel(x(2), y(2), m, anchor);
end