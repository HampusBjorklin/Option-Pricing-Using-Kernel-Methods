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

function xc = ClusterPts(x_min,x_max,N,K,ell,mode)
% Mode can be cheb or sinh
% Assuming that we can have multiple cluster points, and then calling for
% one at a time. Don't really know the effects though.
% If we do Chebyshev, we can let the previouspoint sit in zero, which
% will help. Can we do several periods of Chebyshev?
% Om vi delar in intervallet i först proportionerligt antal punkter och
% sedan fixar varje bit med Cheb eller LeftCheb eller RightCheb.  
%

  i1 = [x_min; K(:)];
  iN = [K(:); x_max];
  do = ones(length(i1),2);
  do(1,1) = 0; do(end,2) = 0;
  if (K(1)==x_min)
    i1 = i1(2:end);
    iN = iN(2:end);
    do = do(2:end,:);
  end
  if (K(end)==x_max)
    i1 = i1(1:end-1);
    iN = iN(1:end-1);
    do = do(1:end-1,:);
  end
  nl = length(i1);

  if (strcmp(mode,'cheb'))
    %
    % Now determine the number of points per interval proportional to length.
    % I think this is OK.
    %
    L = iN-i1
    l = L/sum(L) % Rounding upwards because the endpoints of the interval
    n = ceil(N*l);               % will be duplicated and removed.
    n = max(n,2);
    %
    % For each interval cluster points to one or both endpoints.
    %
    xc=zeros(1,0);


    for k=1:nl
      if sum(do(k,:))==2
        x = ChebPts([i1(k) iN(k)],n(k),ell);
      elseif do(k,1)==1
        x = ChebPtsLeft([i1(k) iN(k)],n(k),ell);
      elseif do(k,2)==1  
        x = ChebPtsLeft(-[iN(k) i1(k)],n(k),ell);
        x = -x(end:-1:1);
      end
      xc = [xc; x(1:end-1)];
    end
    xc = [xc;x(end)];
  else % mode=sinh
    % Now we want to have the cluster points in the middle or one end
    % instead, because I do not know how to make two with this approach.
    % x---y---x---y---x---y---x
    midpts = 0.5*(i1+iN);
    i1 = [x_min; midpts(:)];
    iN = [midpts(:); x_max];
    if (K(1)>i1(2));
      i1 = i1([1 3:end]);
      iN = iN([2:end]);
    end
    if (K(end)<i1(end))
      i1 = i1([1:end-1]);
      iN = iN([1:end-2 end]);
    end
    nl = length(i1);
    L = iN-i1;
    l = L/sum(L); % Rounding upwards because the endpoints of the interval
    n = ceil(N*l);               % will be duplicated and removed.
    n = max(n,2);
    %
    % For each interval cluster to one point
    %
    xc=zeros(1,0);
    for k=1:nl
      x = SinhPts(i1(k),iN(k),n(k),K(k),ell);
      xc = [xc; x(1:end-1)];
    end
    xc = [xc;x(end)];
  end


  
