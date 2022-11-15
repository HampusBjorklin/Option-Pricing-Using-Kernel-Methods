function K = MQ(x,y)
r = norm(x - y, 2);
K = sqrt(1 + eps^2 * r^2);
end

