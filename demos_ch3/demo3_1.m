% Bayesian data analysis
% Aki Vehtari <Aki.Vehtari@aalto.fi>

% Visualise the joint density and marginal densities of the posterior
% of normal distribution with unknown mean and variance.
%
% Note that here pdf's are used for the illustration.  In real
% applications it is better do computations using log-pdfs to maintain
% the numerical stability when magnitude of the density is very small
% or very large.

% data
y=[93 112 122 135 122 150 118 90 124 114];
% sufficient statistics
n=length(y);
s2=var(y);
my=mean(y);

% Factorize the joint posterior p(mu,sigma2|y) to p(sigma2|y)p(mu|sigma2,y)
% Sample from the joint posterior using this factorization

% sample from p(sigma2|y) (sigma2  is a vector 1000 x 1)
sigma2=sinvchi2rand(n-1,s2,1000,1);
% sample from p(mu|sigma2,y) 
mu=my+sqrt(sigma2/n).*randn(size(sigma2));
% display sigma instead of sigma2
sigma=sqrt(sigma2);

clf
% default settings for the figure
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'Position',[200 60 820 640]);
set(gca,'Position',[0.05 0.5 0.5 0.48],'Box','on')

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
fprintf('.')
pause

% plot samples from the joint posterior
% use hold on to plot over the old plot
hold on
plot(mu,sigma,'.','Color',[0 0.5 0]);
hold off
h=legend('Exact contour plot','Samples from the joint posterior',2);
set(h,'FontSize',14)
xlim(tl1)
ylim(tl2)
xlabel('mu')
ylabel('sigma')
% pause
pause

% illustrate the projection to mu
hl=line([mu mu]',[sigma repmat(tl2(1),size(sigma))]','Color','k','LineStyle',':','LineWidth',1);
h=legend('Exact contour plot','Samples from the joint posterior','Projection to mu',2);
set(h,'FontSize',14)
% draw now, then compute something before pause. This way after
% the pause we do not need to wait for those computations
drawnow

% compute the exact marginal density
% multiplication by 1./sqrt(s2/n) is due to the transformation of
% variable z=(x-mean(y))/sqrt(s2/n), see Gelman et al p. 24
pm=tpdf((t1-mean(y))/sqrt(s2/n),n-1)./sqrt(s2/n);
% estimate the marginal density using samples and an ad hoc Gaussian kernel approximation
% kernelp is available from the course web page
pk=kernelp(mu,t1);
% kernelp takes some time, so notify when ready
fprintf('.')
% pause
pause
% change the legend in the first subplot
h=legend('Exact contour plot','Samples from the joint posterior',2);
set(h,'FontSize',14)
% plot the marginal posterior p(mu|y)
axes('Position',[0.05 0.07 0.5 0.37],'Box','on')
plot(t1,pm,t1,pk,'--')
set(gca,'YTick',[])
legend('Exact','Empirical',2)
xlabel('Marginal of mu')
% delete the projection lines
delete(hl)
pause

% illustrate the projection to sigma
set(gcf,'CurrentAxes',ha)
hl=line([mu repmat(tl1(2),size(mu))]',[sigma sigma]','Color','k','LineStyle',':','LineWidth',1);
h=legend('Exact contour plot','Samples from the joint posterior','Projection to sigma',2);
set(h,'FontSize',14)
drawnow

% compute the exact marginal density
% the multiplication by 2*t2 is due to the transformation of
% variable z=t2^2, see Gelman et al p. 24
pm=sinvchi2pdf(t2.^2,n-1,s2)*2.*t2;
% estimate the marginal density using samples and an ad hoc Gaussian kernel approximation
pk=kernelp(sigma,t2);
% kernelp takes some time, so notify when ready
fprintf('.')
% pause
pause
% change the legend in the first subplot
h=legend('Exact contour plot','Samples from the joint posterior',2);
set(h,'FontSize',14)
% plot the marginal distribution p(sigma2|y)
axes('Position',[0.60 0.5 0.37 0.48],'Box','on')
plot(pm,t2,pk,t2,'--')
set(gca,'XTick',[],'YLim',tl2)
legend('Exact','Empirical')
xlabel('Marginal of sigma')
% delete the projection lines
delete(hl)
