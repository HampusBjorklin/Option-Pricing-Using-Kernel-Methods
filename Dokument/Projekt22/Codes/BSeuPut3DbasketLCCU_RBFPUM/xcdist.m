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






