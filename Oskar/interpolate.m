function sf = interpolate(alpha, Repfunc, points, m)
%Interpolate the values of points given alpha and repuducuing kernal
sz = length(points(:,1));
sf = zeros(sz,sz);

for k1 = 1:sz
    for k2 = 1:sz
        for i = 1:N
            sf(k1,k2)= sf(k1,k2) + alpha(i) * (Repfunc([points(k1,1), points(k2,1)], X(i,:), m));
        end
    end
end
end

