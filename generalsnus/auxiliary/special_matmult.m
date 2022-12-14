function [Di] = special_matmult(coeff, mat)
%Multiplies a array of matrices [A, B, C...] with array
% of coefficientns [a, b, c...] s.t the result is 
% a*A + b*B + c*C...
szcoeff = length(coeff);
sz = size(mat);
mat_vec = reshape(mat, numel(mat)/szcoeff, szcoeff);
di = mat_vec*coeff';
Di = reshape(di, sz(1),sz(2)/szcoeff);
end

