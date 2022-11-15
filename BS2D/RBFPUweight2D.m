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

function [phi,phix,phiy,phixx,phixy,phiyx,phiyy] = RBFPUweight2D(x,C,R)
% Compute Wendland functions and its derivatives
% 
% ---------   input   -----------
% x: grid nodes
% C: center nodes
% R: radius of partitions
%
% --------   output  ---------
% phi: Wendland finction
% phix: 1st derivative of Wendland finction in x
% phiy: 1st derivative of Wendland finction in y
% phixx: 2nd derivatiove of Wendland finction in x
% phixy: 2nd mixed derivatiove of Wendland finction
% phiyy: 2nd derivatiove of Wendland finction in y
%
%  (C) Victor Shcherbakov 2014


[phi,phix,phiy,phixx,phixy,phiyx,phiyy] = deal(zeros(size(x,1),1));

ep = 1/R;
dx = x(:,1)-C(1);
dy = x(:,2)-C(2);
r = sqrt(dx.^2 + dy.^2);


phi = (4*ep*r+1).*(1-ep*r).^4;
phix = -20*ep^2*dx.*(1-ep*r).^3;
phiy = -20*ep^2*dy.*(1-ep*r).^3;
phixx = -20*ep*(1-ep*r).^2 + 60*ep^5*dx.^2.*r - 120*ep^4*dx.^2 + 60*ep^3*dx.^2./r;
phixy = 60*ep^5*dx.*dy.*r - 120*ep^4*dx.*dy + 60*ep^3*dx.*dy./r;
phiyx = 60*ep^5*dy.*dx.*r - 120*ep^4*dy.*dx + 60*ep^3*dy.*dx./r;
phiyy = -20*ep*(1-ep*r).^2 + 60*ep^5*dy.^2.*r - 120*ep^4*dy.^2 + 60*ep^3*dy.^2./r;
    
end

















