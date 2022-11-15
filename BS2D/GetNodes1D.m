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

function [xc]=GetNodes1D(x_min,x_max,type,N,K,ell,mode,scale)
%
  switch lower(type)
    case 'uni'
      %
      % Uniform grid with NxN points
      %
      xc = linspace(x_min,x_max,N)';
    case 'cheb'
      xc = ChebPts([x_min,x_max],N);
    case 'cluster'
      xc = scale*ClusterPts(x_min/scale,x_max/scale,N,K/scale,ell,mode);
    otherwise
      error('Requested node type not implemented')
  end
  