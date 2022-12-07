function [X] = getXVector(anchor, n)
%GETXVECTOR Genererar X-vektorn
    dim = length(anchor);
    if prod(size(anchor) == [1, dim])
        anchor = anchor';
    end
    startpoints = anchor.*(ones(dim));
    endpoints = anchor.*(ones(dim));

    startpoints(1,1) = 0;
    endpoints(1,1) = 1;
    for d=2:dim
        startpoints(d, d) = -0.5;
        endpoints(d, d) = 0.5;
    end
    sz = size(startpoints);

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

