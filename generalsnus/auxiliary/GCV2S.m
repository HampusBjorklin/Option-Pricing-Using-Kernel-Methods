function [C] = GCV2S(dim)
% Generates a transformation matrix from S to V
% Done with Gram-Schimt following a vector [1,1,1,1] and easy basis.
% 
A = eye(dim);
A(:,1) = ones(1,dim);

B = zeros(dim,dim);
C = zeros(dim,dim);

for j = 1:dim
    v = A(:,j);
    for i = 1:j-1
        B(i,j) = C(:,i)'*A(:,j);
        v = v -B(i,j)*C(:,i);
    end
B(j,j) = norm(v);
C(:,j) = v/B(j,j);
end
end