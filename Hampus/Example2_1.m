% SYSTEMMATRIX2D 1
%
% Output 3
% B t h e N x N MQ RBF sy s tem ma tr ix
% In pu t 5
% xc N d i s t i n c t c e n t e r s
function B = systemMatrixMQ2d( xc , shape )
[M,N] = size(xc) ;
if N>M, xc = xc' ; end
x = xc (:,1); 
y = xc (:,2);
o = ones(1 ,length(x)) ;
r = sqrt((x*o - ( x*o )').^2 + + (y*o - (y*o )').^2);
B = mq(r, shape);

