function A=RBFmat(phi,ep,r,nprime,dim);

if (nargin==5)  
  A = feval(phi,ep,r,nprime,dim); 
elseif (nargin==4)
  A = feval(phi,ep,r,nprime);   
else
  error('Wrong number of arguments to RBFmat')
end




