function y = sinvchi2pdf(x,nu,s2)
%SINVCHI2_PDF Inverse-Gamma probability density function.
%   Y = SINVCHI2PDF(X,NU,S2) returns the scaled inverse-chi2 probability
%   density function with parameters NU and S2, at the values in X.
%
%   Note: Parametrization as in (Gelman et al, 1996).
%      NU is degrees of freedom
%      S2 is scale

% Copyright (c) 2003 Aki Vehtari

if nargin < 2, 
   error('Requires at least two input arguments.'); 
end

y = log(nu/2).*(nu/2) -gammaln(nu/2) + log(s2)/2*nu - log(x).*(nu/2+1) -nu.*s2/2./x;
y=exp(y);
