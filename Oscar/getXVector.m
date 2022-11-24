function [X] = getXVector(anchor, n, theta)
%GETXVECTOR Roterar korset theta radianer runt ankaret
    startpoints = rotatePoint([[anchor(1); 0], [1; anchor(2)]], anchor, theta);
    endpoints = rotatePoint([[anchor(1); 1], [0; anchor(2)]], anchor, theta);

    sz = size(startpoints);
    dim = 2;
    X = zeros(sz(2)*n+1, dim);
    for i = 1:sz(2)
        for d = 1:dim
            if startpoints(d, i) == endpoints(d, i)
                xi = ones(1, n)*startpoints(d, i);
            else
                pts = ceil(abs((n+1)*(anchor(d)-startpoints(d, i)))...
                    /abs((endpoints(d, i)-startpoints(d, i))));
                l1 = linspace(startpoints(d, i), anchor(d), pts);
                l2 = linspace(anchor(d), endpoints(d, i), n - pts+2);
                xi = [l1(1:end-1) l2(2:end)];
            end
            X(1+(i-1)*n:i*n, d) = xi;
        end
    end
    X(end, :) = anchor;
end

