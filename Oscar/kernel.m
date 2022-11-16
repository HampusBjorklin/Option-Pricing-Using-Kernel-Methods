function [k] = kernel(x, y, m, anchor)
    k = 0;
    if min(x, y) > max(anchor)
        for r = 1:m-1
            k = k + ((x-anchor(1)).^r * (y-anchor(2)).^r)* 1/factorial(r)^2;
        end
        fun = @(z) (x - z).^(m-1) .* 1/factorial(m-1) ...
              .* (y - z).^(m-1) .* 1/factorial(m-1);
        int = integral(fun, max(anchor), min(x, y));
        k = k + int;
    elseif max(x, y) < min(anchor)
        for r = 1:m-1
            k = k + ((anchor(1)-x).^r * (anchor(2)-y).^r)* 1/factorial(r)^2;
        end
        fun = @(z) (z-x).^(m-1) .* 1/factorial(m-1) ...
              .* (z-y).^(m-1) .* 1/factorial(m-1);
        int = integral(fun, max(x, y), min(anchor));
        k = k + int;
        
    else
    end
%     if max(x, y) > anchor && min(x, y) < anchor
%         return
%     end
%     
%     for r = 1:m-1
%         k = k + (max(x-anchor).^r * max(y-anchor).^r)* 1/factorial(r)^2;
%     end
%     fun = @(z) max(x - z, 0).^(m-1) .* 1/factorial(m-1) ...
%           .* max(y - z, 0).^(m-1) .* 1/factorial(m-1);
%     if max(x, y) > anchor 
%         int = integral(fun, anchor, min(x, y));
%     elseif min(x, y) < anchor
%         int = integral(fun, max(x, y), anchor);
%     else
%         disp("test")
%         int = 0;
%     end
%     k = k + int;
end