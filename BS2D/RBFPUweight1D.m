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

function [phi,phix,phixx] = RBFPUweight1D(x,C,R)
% Compute Wendland functions and its derivatives
% 
% ---------   input   -----------
% x: grid nodes
% C: center nodes
% R: radius of partitions
%
% --------   output  ---------
% phi: Wendland finction
% phix: 1st derivative of Wendland finction
% phixx: 2nd derivatiove of Wendland finction
%
%  (C) Victor Shcherbakov 2014

[phi,phix,phixx] = deal(zeros(size(x,1),1));

ep = 1/R;
r = sqrt((x-C).^2)*ep;

phi = (4*r+1).*(1-r).^4;
phix= 4*(r-1).^4 + 4*(4*r+1).*(r-1).^3;
phixx = 32*(r-1).^3 + 12*(4*r+1).*(r-1).^2;
    
end

















