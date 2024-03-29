function rep_kernel = GeneralRepKernelSecondDer(x,y, eps, der_dim, max_order)
%d = dimension
%max_order = maximum allowed order of subsets of the set of dimensions
%x and y are d-dimensional points
%eps is the shape parameter
%der_dim is the dimension we take the derivate in

d = length(x);

if(nargin < 5)
    if d == 2
        max_order = 2;
    else
        max_order = 3;
    end
end

%Function to determine coefficient in front of each derivative 
derivative_coeff = (eps^2)/((eps^2*(x(der_dim)-y(der_dim))^2+1)^(3/2));

%Multiquadric reproducing kernel function
multi = @(a,b) sqrt(1+eps^2*(a-b)^2);

%Cell array of all subset combinations
%s = cell(nr_of_subsets,d);
s = {};

dim = 1:d; %Vector with all dimensions except for the one we look at
dim=dim(dim~=der_dim);

for i=1:max_order-1
    subsets = nchoosek(dim,i);
    for j = 1:length(subsets)
        s(end+1) = {subsets(j,:)};
    end
end

rep_kernel = 1; %+1 is always included

for i=1:length(s)
    arr = cell2mat(s(i));
    arr_len = length(arr);
    coeff = 1;
    for j=1:arr_len
        coord_x = x(arr(j));
        coord_y = y(arr(j));
        coeff = coeff*multi(coord_x,coord_y);

    end
    rep_kernel = rep_kernel + coeff;
end 

rep_kernel = derivative_coeff*rep_kernel;
disp(rep_kernel);
end

