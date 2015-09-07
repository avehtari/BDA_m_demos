% Bayesian data analysis
% Aki Vehtari <Aki.Vehtari@aalto.fi>

% Visualise sampling from the posterior predictive distribution.

% data
y=[93 112 122 135 122 150 118 90 124 114];
% sufficient statistics
n=length(y);
s2=var(y);
my=mean(y);

% Factorize the joint posterior p(mu,sigma2|y) to p(sigma2|y)p(mu|sigma2,y)
% Illustrate the predictive distribution and sampling from it

% sample from the marginal p(sigma2|y)
sigma2=(n-1)*s2./chi2rnd(n-1,1000,1);
% sample from the conditional p(mu|sigma2,y)
mu=my+sqrt(sigma2./n).*randn(size(sigma2));
% display sigma instead of sigma2
sigma=sqrt(sigma2);
% sample from the predictive distribution p(ynew|y)
ynew=randn(size(mu)).*sigma+mu;

clf
% default settings for the figure
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'Position',[200 60 820 640]);
set(gca,'Position',[0.05 0.53 0.43 0.42],'Box','on')
ha=gca;

% For mu compute the density in these points
tl1=[90 150];
t1=linspace(tl1(1),tl1(2),1000);
% For sigma compute the density in these points
tl2=[10 60];
t2=linspace(tl2(1),tl2(2),1000);

% plot the contour plot of the exact posterior
% evaluate density in grid
[T1,T2]=meshgrid(t1,t2);
% note that the following is not normalized, but for plotting
% contours it does not matter
Z=sinvchi2pdf(T2.^2,n-1,s2)*2.*T2.*normpdf(T1,my,T2/sqrt(n));
% linspace(1e-5,max(Z(:)),6) is used to give a vector of linearly
% spaced values at which levels contours are drawn
[cc,hc]=contour(T1,T2,Z,linspace(1e-5,max(Z(:)),6),'b');
set(gca,'XLim',tl1,'YLim',tl2)
xlabel('mu')
ylabel('sigma')

line(mu,sigma,'LineStyle','none','Marker','.','Color',[0 0.5 0]);
yl=[10 60];
set(gca,'XLim',[90 150],'YLim',yl)
xlabel('mu')
ylabel('sigma')
% legend
h=legend('Exact contour plot','Samples from the joint posterior',2);
set(h,'FontSize',14)

% plot ynew
axes('Position',[0.54 0.07 0.43 0.42],'Box','on')
xlabel('ynew')
set(gca,'ytick',[])
set(gca,'XLim',[50 185],'Ylim',[0 0.05])
ha2=gca;

% pick sample from the posterior of mu and sigma
i1=1;
set(gcf,'CurrentAxes',ha)
hl1=line(mu(i1),sigma(i1),'LineStyle','none','Marker','.','Color','r','markersize',20);
h=legend('Exact contour plot','Samples from the joint posterior','A sample from the joint posterior',2);
set(h,'FontSize',14)
pause

set(gcf,'CurrentAxes',ha2)
x=linspace(50,185,1000);
% plot the predictive distribution of ynew given mu and sigma2
hl2=line(x,normpdf(x,my,sqrt(sigma2(i1)/n)),'Color','b','LineStyle','-','LineWidth',2);
% plot a sample form the predictive distribution
hl3=line(ynew(i1),0.0001,'LineStyle','none','Marker','.','Color','r','MarkerSize',20);
h=legend('The pred.dist. given the posterior sample','A sample from the predictive distribution',2);
set(h,'FontSize',14)

pause
h=legend('The pred.dist. given the posterior sample','Samples from the predictive distribution',2);
set(h,'FontSize',14)
delete(hl1)
delete(hl2)
set(hl3,'Color',[0 0.5 0],'Markersize',6)
for i1=2:50
  % first sample sigma2 (here we plot sigma instead)
  set(gcf,'CurrentAxes',ha)
  hl1=line(mu(i1),sigma(i1),'LineStyle','none','Marker','.','Color','r','markersize',20);
  pause(0.1)
  set(gcf,'CurrentAxes',ha2)
  % plot the predictive distribution of ynew fiven mu and sigma2
  hl2=line(x,normpdf(x,my,sqrt(sigma2(i1)/n)),'Color','b','LineStyle','-','LineWidth',2);
  pause(0.1)
  % sample mu given sigma2 
  % plot a sample form the predictive distribution
  hl3=line(ynew(i1),0.0001,'LineStyle','none','Marker','.','Color','r','MarkerSize',20);
  pause(0.1)
  delete(hl1)
  delete(hl2)
  set(hl3,'Color',[0 0.5 0],'Markersize',6)
end
set(gcf,'CurrentAxes',ha)
h=legend('Exact contour plot','Samples from the joint posterior',2);
set(h,'FontSize',14)


pause
% plot the exact marginal
set(gcf,'CurrentAxes',ha2)
pm=tpdf((x-mean(y))/sqrt(s2+s2/n),n-1)./sqrt(s2+s2/n);
h=plot(x,pm);
set(gca,'YTick',[],'YLim',[0 0.03],'XLim',[50 185])
% plot all the samples
hl2=line(ynew,rand(size(ynew))*0.0005,'LineStyle','none','Marker','.','Color',[0 0.5 0]);
h=legend([h hl2],'Exact predictive distribution','Samples from the predictive distribution');
set(h,'FontSize',14)
xlabel('ynew')

