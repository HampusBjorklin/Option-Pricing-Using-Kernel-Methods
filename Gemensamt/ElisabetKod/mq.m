 %Copyright (C) 2015 Victor Shcherbakov

    %This file is part of BENCHOP.
    %BENCHOP is free software: you can redistribute it and/or modify
    %it under the terms of the GNU General Public License as published by
    %the Free Software Foundation, either version 3 of the License, or
    %(at your option) any later version.

    %BENCHOP is distributed in the hope that it will be useful,
    %but WITHOUT ANY WARRANTY; without even the implied warranty of
    %MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %GNU General Public License for more details.

    %You should have received a copy of the GNU General Public License
    %along with BENCHOP. If not, see <http://www.gnu.org/licenses/>.

function [phi]=mq(epsil,r,nprime,dim)
%
% NPRIME is a string defining which operator to use on the basis function
%
% DIM is the dimension for the partial derivative if nprime is '0','1'...,'4'
% DIM(1:2) are the dimensions for the mixed second derivative if nprime
% is 'm2'
% and DIM is the number of space dimensions if nprime is 'L' or 'L2'
%
%
% Assume a one-dimensional problem if no dimension is given
% 
if nargin<=3
  dim=1;
end
%
% For the case of L or L2 operators, we need to know the number of dimensions
%
if (nprime(1)=='L')
  if (size(r,3)==1)
    nd=dim;
  else
    nd=size(r,3)-1;
  end
end
%
% For the mixed derivative, the dimensions must be given even in 2D
%
if (nprime(1)=='m')
  if (length(dim)~=2)
    error('For the mixed derivative, dim=dim(1:2)')
  elseif (dim(1)==dim(2))
    error('For mixed derivatives, dim(1) must be other than dim(2)')
  end  
end
%
% epsil can be either just one value or a vector of N values
%
esz = size(epsil);
if (prod(esz)~=1)
  if (min(esz)==1 & max(esz)==size(r,2))
    %
    % Make epsil into a matrix with constant columns
    %
    epsil = ones(size(r,1),1)*epsil(:).';
  else
    error('The size of epsil does not match the columns in r')
  end
end

phi=zeros(size(r,1),size(r,2));

tmp = 1+(epsil.*sq(r(:,:,1))).^2;

if nprime(1)=='0'
  phi = sqrt(tmp);

elseif nprime(1)=='1'
  phi = epsil.^2.*sq(r(:,:,dim+1))./sqrt(tmp);

elseif nprime(1)=='2'
  phi = epsil.^2./sqrt(tmp) - epsil.^4.*sq(r(:,:,dim+1)).^2.*tmp.^-1.5;

elseif nprime(1)=='3'
  phi = -3*epsil.^4.*sq(r(:,:,dim+1))   .*tmp.^-1.5 + ...
         3*epsil.^6.*sq(r(:,:,dim+1)).^3.*tmp.^-2.5;

elseif nprime(1)=='4'
  phi = -3*epsil.^4.*tmp.^-1.5 +18*epsil.^6.*sq(r(:,:,dim+1)).^2.*tmp.^-2.5 - ...
                              15*epsil.^8.*sq(r(:,:,dim+1)).^4.*tmp.^-3.5;

elseif nprime(1)=='L' & length(nprime)==1
  phi = nd*epsil.^2.*tmp.^-0.5 - epsil.^4.*sq(r(:,:,1)).^2.*tmp.^-1.5;

elseif nprime(1:2)=='L2'
  phi = -nd*(nd+2)*epsil.^4.*tmp.^-1.5 + ...
         6*(nd+2)*epsil.^6.*sq(r(:,:,1)).^2.*tmp.^-2.5 - ...
         15*epsil.^8.*sq(r(:,:,1)).^4.*tmp.^-3.5; 

elseif nprime(1:2)=='m2'
  phi = -epsil.^4.*sq(r(:,:,dim(1)+1)).*sq(r(:,:,dim(2)+1)).*tmp.^-1.5; 
  
else
  error('Error in input argument nprime to function mq')
end 

function r=sq(r)
 r=squeeze(r);
