function [rotatedPoint] = rotatePoint(p, c, theta)
%ROTATEPOINT Summary of this function goes here
%   Detailed explanation goes here
r = @(p, c, theta) [cos(theta), sin(theta); -sin(theta), cos(theta)]*(p-c)+c;
rotatedPoint = r(p, c, theta);
sz = size(p);

for i = 1:sz(2)
    while sum(rotatedPoint(:, i)>[1; 1]) > 0 || sum(rotatedPoint(:, i)<[0; 0]) > 0 % Outside the square (0, 0) - (1, 1)
        rotatedPoint(:, i) = rotatedPoint(:, i) + 0.01*(c-rotatedPoint(:, i));
    end
end

