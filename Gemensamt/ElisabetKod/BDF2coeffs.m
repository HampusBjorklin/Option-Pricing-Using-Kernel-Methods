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

function [k,beta0,beta1,beta2]=BDF2coeffs(T,M)
%
% T is the final time, M is the number of timesteps.  
% beta0 is constant throughout the time-stepping.
% beta1 and beta2 are vectors with coefficients for the M steps
%  
  c(1)=1;
  k(1)=1; 
  for step=2:M
    omega(step,1)=c(step-1) - 0.5 + 0.5*sqrt(4*c(step-1)^2+1);
    k(step,1) = omega(step)*k(step-1);
    c(step,1) = (1+omega(step))/(1+2*omega(step));
  end
  %
  % Scale the steps to fit the final time
  %
  k = T/sum(k)*k;
  %
  % Compute the coefficients for the method. (Initial values are correct.)
  %
  beta0=k(1);
  beta1 = (1+omega).^2./(1+2*omega);
  beta2 =     omega.^2./(1+2*omega);
