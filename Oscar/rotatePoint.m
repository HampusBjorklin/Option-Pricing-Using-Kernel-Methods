function [rotatedPoint] = rotatePoint(p, c, theta)
%ROTATEPOINT Summary of this function goes here
%   Detailed explanation goes here
r = @(p, c, theta) [cos(theta), sin(theta); -sin(theta), cos(theta)]*(p-c)+c;
rotatedPoint = r(p, c, theta);

okDist = 0.025; % point is ok if it is this close to the boundary
sz = size(p);
for i = 1:sz(2)
    if sum(rotatedPoint(:, i)>[1; 1]) > 0 || sum(rotatedPoint(:, i)<[0; 0]) > 0
        while sum(rotatedPoint(:, i)>[1; 1]) > 0 || sum(rotatedPoint(:, i)<[0; 0]) > 0 % Outside the square (0, 0) - (1, 1)
            rotatedPoint(:, i) = rotatedPoint(:, i) + 0.01*(c-rotatedPoint(:, i));
        end
    elseif sum(rotatedPoint(:, i)<[1-okDist; 1-okDist]) == 2 && sum(rotatedPoint(:, i)>[okDist; okDist]) == 2
        while sum(rotatedPoint(:, i)<[1-okDist; 1-okDist]) == 2 && sum(rotatedPoint(:, i)>[okDist; okDist]) == 2
            rotatedPoint(:, i) = rotatedPoint(:, i) - 0.05*(c-rotatedPoint(:, i));
        end
    end
end

