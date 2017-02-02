% Bayesian data analysis
% Aki Vehtari <Aki.Vehtari@aalto.fi>

% Visualise the marginal distribution of mu as a mixture of normals.

% data
y=[93 112 122 135 122 150 118 90 124 114];
% sufficient statistics
n=length(y);
s2=var(y);
my=mean(y);

% Factorize the joint posterior p(mu,sigma2|y) to p(sigma2|y)p(mu|sigma2,y)
% Illustrate marginal of mu as a scale mixture of normals, where
% mixing parameter is distributed as p(sigma2|y)

% sample from p(sigma2|y)
sigma2=(n-1)*s2./chi2rnd(n-1,1000,1);
% sample from p(mu|sigma2,y)
mu=my+sqrt(sigma2/n).*randn(size(sigma2));
% display sigma instead of sigma2
sigma=sqrt(sigma2);

clf
% default settings for the figure
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'Position',[200 60 820 640]);

% For mu compute the density in these points
tl1=[90 150];
t1=linspace(tl1(1),tl1(2),1000);
% For sigma compute the density in these points
tl2=[10 60];
t2=linspace(tl2(1),tl2(2),1000);

% Create axes for the joint distribution
set(gca,'Position',[0.05 0.5 0.5 0.48],'Box','on')
% store axes handle
ha=gca;
% plot the contour plot of the exact posterior
% evaluate density in grid
[T1,T2]=meshgrid(t1,t2);
% note that the following is not normalized, but for plotting
% the contours it does not matter
Z=sinvchi2pdf(T2.^2,n-1,s2)*2.*T2.*normpdf(T1,my,T2/sqrt(n));
% linspace(1e-5,max(Z(:)),6) is used to give a vector of linearly
% spaced values at which levels contours are drawn
[cc,hc]=contour(T1,T2,Z,linspace(1e-5,max(Z(:)),6),'b');
set(gca,'XLim',tl1,'YLim',tl2)
xlabel('mu')
ylabel('sigma')
% legend
h=legend('Exact contour plot','location','northwest');
set(h,'FontSize',14)

% compute the exact marginal density
% the multiplication by t2 is due to the transformation of
% variable z=t2^2, see BDA3 p. 21
pm=sinvchi2pdf(t2.^2,n-1,s2)*2.*t2;
% plot the marginal distribution of sigma
axes('Position',[0.60 0.5 0.37 0.48],'Box','on')
plot(pm,t2)
set(gca,'XTick',[],'YLim',tl2)
xlabel('Marginal of sigma')

% Create axes for the marginal of mu
ha2=axes('Position',[0.05 0.07 0.5 0.37],'Box','on');
xlabel('Marginal of mu')

fprintf('.')
pause

% illustrate sampling
x=linspace(tl1(1),tl2(2),1000);
for i1=1:50
  % first sample sigma2 (here we plot sigma instead)
  set(gcf,'CurrentAxes',ha)
  hl1=line(get(gca,'XLim'),[sigma(i1) sigma(i1)]','Color','k','LineStyle','--','LineWidth',2);
  pause(0.1)
  % plot conditional distribution of mu given sigma2
  % (scaling of pdf is just for illustration)
  hl2=line(t1,sigma(i1)+normpdf(t1,my,sqrt(sigma2(i1)/n))*100,'Color',[0 0.5 0],'LineStyle','--','LineWidth',2);
  if i1==1
    % plot legend only once
    h=legend([hc, hl1, hl2],'Exact contour plot','Sample from the marginal of sigma','Conditional distribution of mu','location','northwest');
    set(h,'FontSize',14)
  end
  set(gcf,'CurrentAxes',ha2)
  hl3(i1)=line(t1,normpdf(t1,my,sqrt(sigma2(i1)/n)),'Color',[0 0.5 0],'LineStyle','--','LineWidth',1);
  set(gca,'YTick',[],'YLim',[0 0.1])
  if i1==1
    % plot legend only once
    h=legend('Samples of conditionals of mu','location','northwest');
    set(h,'FontSize',14)
  end
  pause(0.1)
  % sample mu given sigma2 
  set(gcf,'CurrentAxes',ha)
  hs=line(mu(i1),sigma(i1),'LineStyle','none','Marker','.','Color',[0 .5 0]);
  if i1==1
    % plot legend only once
    h=legend([hc, hl1, hl2, hs],'Exact contour plot','Sample from the marginal of sigma','Conditional distribution of mu','Samples from the joint posterior','location','northwest');
    set(h,'FontSize',14)
  end
  pause(0.1)
  delete(hl1)
  delete(hl2)
end
h=legend([hc, hs],'Exact contour plot','Samples from the joint posterior','location','northwest');
set(h,'FontSize',14)
    
% compute the marginal distribution using mixture of conditional distributions
% in this case conditional distributions are normal distributions
for i1=1:500
  pps(i1,:)=normpdf(t1,my,sqrt(sigma2(i1)/n));
end
pp=mean(pps);
fprintf('.')
pause
% plot all 1000 samples
hs=line(mu,sigma,'LineStyle','none','Marker','.','Color',[0 .5 0]);
set(gca,'XLim',tl1,'YLim',tl2)
% change axes
set(gcf,'CurrentAxes',ha2)
% compute the exact marginal (t-distribution)
pm=tpdf((t1-mean(y))/sqrt(s2/n),n-1)./sqrt(s2/n);
% plot the exact marginal and average of sampled conditionals
% note that average of sampled conditionals is mixture
% of normal distrbutions
plot(t1,pm,t1,pp)
set(gca,'YTick',[],'YLim',[0 0.1])
% legend
h=legend('Exact marginal','Average of sampled conditionals','location','northwest');
set(h,'FontSize',14);
