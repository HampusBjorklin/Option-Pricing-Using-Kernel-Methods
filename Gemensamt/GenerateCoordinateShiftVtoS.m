function [A] = GenerateCoordinateShiftVtoS(dim)
%GenerateCoordinateShiftVtoS Does the inverse of desired coordinate shift in dimension dim
%   Shifts the interesting diagonal, expressed as the first coordinate in V, 
%   back to S

A = eye(dim);
for d = 2:dim
    Aloc = eye(dim);
    Aloc(1,1) = 1/2;
    Aloc(d,1) = 1/2;
    Aloc(1, d) = 1/2;
    Aloc(d,d) = -1/2;
    A = A * Aloc;
end

end
