function [k] = kernal(x, y, m)
   k = 0;
   for r = 1:m-1
       k = k + (x.^r * y.^r)* 1/factorial(r)^2;
   end
   fun = @(z) max(x - z, 0).^(m-1) .* 1/factorial(m-1) ...
          .* max(y - z, 0).^(m-1) .* 1/factorial(m-1);
   int = integral(fun, 0, 1);
   k = k + int;
end