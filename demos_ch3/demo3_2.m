% Bayesian data analysis
% Aki Vehtari <Aki.Vehtari@aalto.fi>

% Visualise factored sampling and the corresponding marginal and
% conditional density.

% data
y=[93 112 122 135 122 150 118 90 124 114];
% sufficient statistics
n=length(y);
s2=var(y);
my=mean(y);

% Factorize the joint posterior p(mu,sigma2|y) to p(sigma2|y)p(mu|sigma2,y)
% Sample from the joint posterior using this factorization

% sample from p(sigma2|y) (sigma2  is a vector 1000 x 1)
sigma2=(n-1)*s2./chi2rnd(n-1,1000,1);
% sample from p(mu|sigma2,y)
mu=my+sqrt(sigma2./n).*randn(size(sigma2));
% display sigma instead of sigma2
sigma=sqrt(sigma2);

clf
% default settings for the figure
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'Position',[200 60 820 640]);

% Create axes for the joint distribution
set(gca,'Position',[0.05 0.5 0.5 0.48],'Box','on')
ha=gca;

% For mu compute the density in these points
tl1=[90 150];
t1=linspace(tl1(1),tl1(2),1000);
% For sigma compute the density in these points
tl2=[10 60];
t2=linspace(tl2(1),tl2(2),1000);

% plot the contour plot of the exact posterior
% evaluate density in a grid
[T1,T2]=meshgrid(t1,t2);
% note that the following is not normalized, but for plotting
% contours it does not matter
Z=sinvchi2pdf(T2.^2,n-1,s2)*2.*T2.*normpdf(T1,my,T2/sqrt(n));
% linspace(1e-5,max(Z(:)),6) is used to give a vector of linearly
% spaced values at which levels contours are drawn
contour(T1,T2,Z,linspace(1e-5,max(Z(:)),6),'b');
% store axes handle
ha=gca;
% legend
h=legend('Exact contour plot',2);
set(h,'FontSize',14)
%
set(gca,'XLim',tl1,'YLim',tl2)
xlabel('mu')
ylabel('sigma')

% compute the exact marginal density of sigma
% the multiplication by t1 is due to the transformation of
% variable z=t2^2, see BDA3 p. 21
pm=sinvchi2pdf(t2.^2,n-1,s2)*2.*t2;
% plot the marginal distribution of sigma
axes('Position',[0.60 0.5 0.37 0.48],'Box','on')
plot(pm,t2)
set(gca,'XTick',[],'YLim',tl2)
xlabel('marginal of sigma')
fprintf('.')
pause

% illustrate the sampling
% note that actual sampling was done above, and below we just
% illustrate using those samples
set(gcf,'CurrentAxes',ha)
xlabel('mu')
ylabel('sigma')
x=linspace(tl1(1),tl1(2),1000);
% first illustrate step by step one sample
i1=1;
% first sample sigma2 (here we plot sigma instead)
hl1=line(tl1,[sigma(i1) sigma(i1)]','Color','k','LineStyle','--','LineWidth',2);
h=legend('Exact contour plot','Sample from the marginal of sigma',2);
set(h,'FontSize',14)
pause
% plot the conditional distribution of mu given sigma2
% (scaling of pdf is just for illustration)
hl2=line(x,sigma(i1)+normpdf(x,my,sqrt(sigma2(i1)/n))*100,'Color',[0 0.5 0],'LineStyle','--','LineWidth',2);
h=legend('Exact contour plot','Sample from the marginal of sigma','Conditional distribution of mu',2);
set(h,'FontSize',14)
pause
% sample mu given sigma2 
hs=line(mu(i1),sigma(i1),'LineStyle','none','Marker','.','Color',[0 0.5 0],'MarkerSize',30);
h=legend('Exact contour plot','Sample from the marginal of sigma','Conditional distribution of mu','Sample from the joint posterior',2);
set(h,'FontSize',14)
pause
set(hs,'MarkerSize',6)
h=legend('Exact contour plot','Sample from the marginal of sigma','Conditional distribution of mu','Samples from the joint posterior',2);
set(h,'FontSize',14)
% delete the lines showing sigma sample and conditional of mu
delete(hl1)
delete(hl2)
for i1=2:100
  % first sample sigma2 (here we plot sigma instead)
  hl1=line(xlim,[sigma(i1) sigma(i1)]','Color','k','LineStyle','--','LineWidth',2);
  pause(0.1)
  % plot the conditional distribution of mu given sigma2
  % (scaling of pdf is just for illustration)
  hl2=line(x,sigma(i1)+normpdf(x,my,sqrt(sigma2(i1)/n))*100,'Color',[0 0.5 0],'LineStyle','--','LineWidth',2);
  pause(0.1)
  % sample mu given sigma2 
  line(mu(i1),sigma(i1),'LineStyle','none','Marker','.','Color',[0 0.5 0])
  pause(0.1)
  % delete lines showing sigma sample and conditional of mu
  delete(hl1)
  delete(hl2)
end
h=legend('Exact contour plot','Samples from the joint posterior',2);
set(h,'FontSize',14)
