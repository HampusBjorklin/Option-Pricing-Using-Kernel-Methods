function rep_kernel = GeneralRepKernelSecondMixed(x,y, eps, der_dim, max_order)
%d = dimension
%max_order = maximum allowed order of subsets of the set of dimensions
%x and y are d-dimensional points
%eps is the shape parameter
%der_dim is an array with the dimensions we take the derivate in

d = size(x,2);

if(nargin < 5)
    if d == 2
        max_order = 2;
    else
        max_order = 3;
    end
end

%Function to determine coefficient in front of each derivative 
derivative_coeff = (eps^4)*(x(der_dim(1))-y(der_dim(1)))/(sqrt(eps^2*(x(der_dim(1))-y(der_dim(1)))^2+1)) *...
    (x(der_dim(2))-y(der_dim(2)))/(sqrt(eps^2*(x(der_dim(2))-y(der_dim(2)))^2+1));

%Multiquadric reproducing kernel function
multi = @(a,b) sqrt(1+eps^2*(a-b).^2);

%Cell array of all subset combinations
%s = cell(nr_of_subsets,d);
s = {};

dim = 1:d; %Vector with all dimensions except for the one we look at
dim=dim(dim~=der_dim(1));
dim=dim(dim~=der_dim(2));

for i=1:max_order-2
    subsets = nchoosek(dim,i);
    for j = 1:size(subsets,1)
        s(end+1) = {subsets(j,:)};
    end
end

if isempty(s) && max_order<2
    rep_kernel = zeros(size(x,1), 1);
else
   rep_kernel = ones(size(x,1), 1); %+1 is always included
    
    for i=1:length(s)
        arr = cell2mat(s(i));
        arr_len = length(arr);
        coeff = 1;
        for j=1:arr_len
            coord_x = x(:,arr(j));
            coord_y = y(:,arr(j));
            coeff = coeff.*multi(coord_x,coord_y);
    
        end
        rep_kernel = rep_kernel + coeff;
    end 
    
    rep_kernel = derivative_coeff*rep_kernel;
end

end

