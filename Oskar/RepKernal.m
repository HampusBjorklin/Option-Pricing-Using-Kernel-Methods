function [K] = RepKernal(x, y , m)
% Compute the repuducing kernal for Lambda = {null, {1}, {2}}
% X = [x1,..,xn]; y = [y1,...y]
    K = 1 + kernal(x(1), y(1), m) +  kernal(x(2), y(2), m);
end