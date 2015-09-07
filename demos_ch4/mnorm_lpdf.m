function y = mnorm_lpdf(x,mu,sigma)
%NORMPDF Multivariate-Normal log-probability density function (lpdf).
%   Y = MNORM_LPDF(X,MU,SIGMA) Returns the log of the multivariate-normal
%    pdf with mean, MU, and covariance matrix, SIGMA, at the values in X.
%
%   Default values for MU and SIGMA are 0 and 1 respectively.

if nargin < 3, 
  sigma = 1;
end

if nargin < 2;
  mu = 0;
end

if nargin < 1, 
  error('Requires at least one input argument.');
end

[m,n]=size(x);
if m~=1
  error('Vector x must be 1xN.')
end
x=x-mu;
y=-0.5*(sum(sum(inv(sigma).*(x'*x))) +log(det(sigma)) +log(2*pi)*n);
