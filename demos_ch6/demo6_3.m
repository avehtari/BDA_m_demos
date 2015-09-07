% Bayesian data analysis
% Aki Vehtari <Aki.Vehtari@aalto.fi>

% Posterior predictive checking
% Light speed example with poorly chosen test statistic

clf
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')

y=load('light.txt');
% sufficient statistics
n=length(y);
s2=var(y);
s=sqrt(s2);
my=mean(y);

% Second example of replications
for i=1:1000
  pps(:,i)=trnd(n-1,n,1).*sqrt(1+1/n).*s+my;
end
% Use variance as a test statistic
% This is poor choice since it corresponds directly to
% variance parameter in the model which has been fitted
% to the data.
pp=var(pps);
clf
% Plot the distribution of the test statistics (20 = number of bars in histogram)
hist(pp,20);
% Plot the test statistic for the original data set
h1=get(gca,'Children');
h2=line([var(y) var(y)],ylim);
legend([h2 h1],'Variance of original data','Variances of replicated datasets')
title({'Light speed example with poorly chosen test statistic','Pr(T(yrep,\theta)\leq T(y,\theta)|y)=0.42'})

