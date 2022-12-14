function rep_kernel = GeneralRepKernel(x,y, eps, max_order)
%d = dimension
%max_order = maximum allowed order of subsets of the set of dimensions
%max_order must be at least d-1 large
%x and y are d-dimensional points
%eps is the shape parameter

%Cell array of all subset combinations
%s = cell(nr_of_subsets,d);
d = size(x,2);

if(nargin < 4)
    if d == 2
        max_order = 2;
    else
        max_order = 3;
    end
end

s = {};

dim = 1:d; %Vector with all dimensions
for i=1:max_order
    subsets = nchoosek(dim,i);
    for j = 1:size(subsets,1)
        s(end+1) = {subsets(j,:)};
    end
end

%Multiquadric reproducing kernel function
multi = @(a,b) sqrt(1+eps^2*(a-b).^2);

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

end

