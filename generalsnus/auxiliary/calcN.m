function N = calcN(dim, maxorder, n)
N = 1;
for i = 1:maxorder
   N = N + nchoosek(dim,i)*n^i; 
    
end