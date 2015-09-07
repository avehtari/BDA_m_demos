function [e,g,H] = bioassayfun(w,x,y,n)
a=w(1);b=w(2);
% these help using chain rule in derivation
t=a+b.*x;
et=exp(t);
z=et./(1+et);
% negative log posterior (error function to be minimized)
e=-sum(y.*log(z)+ (n-y).*log(1-z));
% Ex 4.2: change these to test your gradients and Hessian
g=NaN;
H=NaN;
