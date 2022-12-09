function [A] = GenerateCoordinateShift(dim)
%GenerateCoordinateShift Does the desired coordinate shift in dimension dim
%   Shifts the x-axis to be the interesting diagonal, and all the
%   coordinates will be orthogonal

A = eye(dim);
for d = 2:dim
    Aloc = eye(dim);
    Aloc(1,1) = 1;
    Aloc(d,1) = 1;
    Aloc(1, d) = 1;
    Aloc(d,d) = -1;
    A = A * Aloc;
end

end
