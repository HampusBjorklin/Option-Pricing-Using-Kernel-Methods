function [C,P]=BS(S0,K,T,r,mu,sigma);
% function [C,P]=BS(S0,K,T,r,mu,sigma);
%
% ----------input---------------------------
%
%      S0: intial/current stock price
%       K: strike price
%       T: time to maturity
%       r: continously compounded short rate
%   sigma: volatility in Black-Scholes model
%  Note that for S0,K,T and r that they should either be 1-dim
%  or match the dimension of P

% (C) 2006 Magnus Wiktorsson
d1=1./(sigma.*sqrt(T)).*(log(S0./K)+mu.*T)+sigma.*sqrt(T);
d2=d1-sigma.*sqrt(T);
C=exp(T.*(sigma.^2/2+mu-r)).*S0.*normcdf(d1)-exp(-r.*T).*K.*normcdf(d2);
P=C+K.*exp(-r.*T)-S0.*exp(T.*(sigma.^2/2+mu-r));
