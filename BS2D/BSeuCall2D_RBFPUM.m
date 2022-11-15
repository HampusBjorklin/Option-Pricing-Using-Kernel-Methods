 %Copyright (C) 2015 Victor Shcherbakov, modified 2022 Elisabeth Larsson

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

function U = BSeuCall2D_RBFPUM(S,K,T,r,sig1,sig2,rho,N,M,ep,Np,del);

% The code expects a column vector, so if S is in a row, fix this
%
%---------   input  ----------------
% S(1:np,1:2) - points where we want to evaluate
% T - time to maturity
% r - risk free interest rate
% sig1 - volatility on 1st asset
% sig2 - volatility on 2nd asset
% rho - correlation
% N - The number of points in each asset dimensions  
% M - the number of time steps
% ep - the shape parameter
% Np - the number of patches in one dimension
% del - relative patch overlap  
%---------  output  ----------------
% U - option value(s)
%
%  (C) Victor Shcherbakov & Elisabeth Larsson 2014

    %
    % First scale everything to the unit square. This is usually beneficial.
    % For baskets scaling is harmless, but may need checking otherwise
    % Use a domain size of 4K
    %
  dim = 2;
  smax = dim*4*K;
  K0 = K;
  S = S/smax; 
  K = K0/smax;   
  %
  % Make function handles for initial condition and boundary conditions 
  %
  payoff = @(S,K,r,t) max(0.5*(S(:,1)+S(:,2))-K*exp(-r*t),0);
  nearbc = @(S,K,r,t) zeros(size(S,1),1);
  farbc = @(S,K,r,t) payoff(S,K,r,t);
  %
  % Define the distances where the boundary conditions are applied
  %
  neard = 0; % Only one point unless we use the log transform
  fard = 1; % The upper right triangle of the domain.
  % We could actually discard the far points for efficiency, but we keep
  % them for now
  % 
  % Define the discretization
  %
  phi = 'mq';
  x = linspace(0,1,N);
  [xx,yy] = meshgrid(x);
  xc = [xx(:) yy(:)];
  

  U=EuRBFPU2D(S,K,T,r,sig1,sig2,rho,payoff, ...
	      nearbc,neard,farbc,fard,phi,ep,xc,M,Np,Np,del);
  %
  % Scale back, note that the payoff was scaled by 1/smax/K0
  % 
  U = smax*U;


