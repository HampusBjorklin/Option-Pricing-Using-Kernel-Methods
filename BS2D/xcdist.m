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

function r=xcdist(x,c,all);
% 
% Evaluation/collocation points and center points are not always the
% same. We compute r=||xi-cj|| together with componentwise signed differences.
%
if (nargin==1)
  all=0;
  c=x;
elseif (nargin==2)
  all=0;
end
  
[np,nd] = size(x);
nc = size(c,1);
nr = 1 + all*nd;

% r contains r and each ri, i=1...nd
r = zeros(np,nc,nr);
for d=1:nd
  [pi,pj] = meshgrid(c(:,d),x(:,d));
  r(:,:,1) = r(:,:,1) + (pi-pj).^2;
  if (all)
    r(:,:,d+1) = pj-pi;
  end
end
r(:,:,1) = sqrt(r(:,:,1));






