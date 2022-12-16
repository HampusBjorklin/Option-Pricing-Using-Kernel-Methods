function [X] = getXVector(anchor, n, maxorder)
%GETXVECTOR Genererar X-vektorn
    if nargin < 3
        maxorder = 1;
    end
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

    if maxorder == 2
        X2 = ones((n^2 + 2*n) * nchoosek(dim, 2) + 1, dim).*anchor';
        pairs = nchoosek(1:dim, maxorder);
        for p=1:size(pairs, 1)
            i = pairs(p, 1);
            j = pairs(p, 2);
            cor1 = [X(1 + (i-1)*n:i*n, i); X(end, i)];
            cor2 = [X(1 + (j-1)*n:j*n, j); X(end, j)];
            [grid1, grid2] = meshgrid(cor1, cor2);
            grid1 = reshape(grid1, [], 1);
            grid2 = reshape(grid2, [], 1);
            X2(1 + (p-1)*(n^2 + 2*n):p*(n^2 + 2*n), [i,j]) = [grid1(1:end-1), grid2(1:end-1)];
        end
        X = X2;
    elseif maxorder == 3
        X3 = ones((n^3 + 3*n^2 + 3*n) * nchoosek(dim, 2) + 1, dim).*anchor';
        trio = nchoosek(1:dim, maxorder);
        for p=1:size(trio, 1)
            i = trio(p, 1);
            j = trio(p, 2);
            k = trio(p, 3);
            cor1 = [X(1 + (i-1)*n:i*n, i); X(end, i)];
            cor2 = [X(1 + (j-1)*n:j*n, j); X(end, j)];
            cor3 = [X(1 + (k-1)*n:k*n, k); X(end, k)];
            [grid1, grid2, grid3] = meshgrid(cor1, cor2, cor3);
            grid1 = reshape(grid1, [], 1);
            grid2 = reshape(grid2, [], 1);
            grid3 = reshape(grid3, [], 1);
            X3(1 + (p-1)*(n^3 + 3*n^2 + 3*n):p*(n^3 + 3*n^2 + 3*n), [i,j, k]) = [grid1(1:end-1), grid2(1:end-1), grid3(1:end-1)];
        end
        X = X3;
    end

end

