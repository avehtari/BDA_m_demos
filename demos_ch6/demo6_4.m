% Bayesian data analysis
% Aki Vehtari <Aki.Vehtari@aalto.fi>

% Marginal posterior predictive checking
% Light speed example

clf
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')

y=load('light.txt');
% sufficient statistics
n=length(y);
s2=var(y);
s=sqrt(s2);
my=mean(y);

% tail area probabilities of marginal predictive distributions
Ty=t_cdf(y,n-1,my,sqrt(1+1/n).*s);
hist(Ty,0.05:.05:.95)
title({'Light speed example','distribution of marginal posterior p-values'})
% Histogram should look uniformish, if the posterior predictive
% distributions are well calibrated (this test is more reliable with
% cross-validation predictive distributions).

