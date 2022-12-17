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
                pts = ceil(abs((n)*(anchor(d)-startpoints(d, i)))...
                    /abs((endpoints(d, i)-startpoints(d, i))));
                l1 = linspace(startpoints(d, i), anchor(d), pts+1);
                l2 = linspace(anchor(d), endpoints(d, i), n - pts+1);
                xi = [l1(1:end-1) l2(2:end)];
            end
            X(1+(i-1)*n:i*n, d) = xi;
        end
    end
    X(end, :) = anchor;
    if maxorder == 2
        sz = n^2;
        X2 = ones((sz) * nchoosek(dim, maxorder), dim).*anchor';
        pairs = nchoosek(1:dim, maxorder);
        i = 0;
        for pair=pairs'
            i = i +1;
            [xx1, yy1] = ndgrid(X(1 + (pair(1)-1)*n:pair(1)*n, pair(1)),...
                                  X(1 + (pair(2)-1)*n:pair(2)*n, pair(2)));

            X2(1 + (i-1)*sz:(i)*sz,pair) = [xx1(:), yy1(:)];
        end
        X = [X;X2];
    elseif maxorder == 3
        sz = n^3;
        X3 = ones((sz) * nchoosek(dim, maxorder), dim).*anchor';
        X2 = getXVector(anchor,n, maxorder -1);
        pairs = nchoosek(1:dim, maxorder);
        i = 0;
        for pair=pairs'
            i = i +1;
            [xx1, yy1, zz1] = ndgrid(X(1 + (pair(1)-1)*n:pair(1)*n, pair(1)),...
                                  X(1 + (pair(2)-1)*n:pair(2)*n, pair(2)),...
                                  X(1 + (pair(3)-1)*n:pair(3)*n, pair(3)));

            X3(1 + (i-1)*sz:(i)*sz,pair) = [xx1(:), yy1(:), zz1(:)];
        end
        X= [X2;X3];
    elseif maxorder == 4
        sz = n^4;
        X4 = ones((sz) * nchoosek(dim, maxorder), dim).*anchor';
        X3 = getXVector(anchor,n, maxorder -1);
        pairs = nchoosek(1:dim, maxorder);
        i = 0;
        for pair=pairs'
            i = i +1;
            [xx1, yy1, zz1, vv1] = ndgrid(X(1 + (pair(1)-1)*n:pair(1)*n, pair(1)),...
                                  X(1 + (pair(2)-1)*n:pair(2)*n, pair(2)),...
                                  X(1 + (pair(3)-1)*n:pair(3)*n, pair(3)),...
                                  X(1 + (pair(4)-1)*n:pair(4)*n, pair(4)));

            X3(1 + (i-1)*sz:(i)*sz,pair) = [xx1(:), yy1(:), zz1(:), vv1(:)];
        end
        X= [X3;X4];
    end
end

