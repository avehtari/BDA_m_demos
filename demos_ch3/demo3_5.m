% Bayesian data analysis
% Aki Vehtari <Aki.Vehtari@aalto.fi>

% Estimating the speed of light using normal model BDA3 p. 66

% data
y=load('light.txt');
% sufficient statistics
n=length(y);
s2=var(y);
my=mean(y);
% sample from p(sigma2|y)
sigma2=(n-1)*s2./chi2rnd(n-1,1000,1);
% sample from p(mu|sigma2,y)
mu=my+sqrt(sigma2/n).*randn(size(sigma2));
% display sigma instead of sigma2
sigma=sqrt(sigma2);

clf
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'Position',[600 60 700 640]);

% histogram of data
subplot(2,1,1)
set(gca,'box','on')
hist(y,-43:2:41)
set(gca,'YTick',[],'XTick',-40:20:40)
title('Newcomb''s measurements')
xlim([-45 45])

pause

% the posterior distribution of mu
subplot(2,1,2)
set(gca,'box','on')
x=linspace(18.5,33.5,1000);
% the exact posterior distribution (t-distribution)
pm=tpdf((x-mean(y))/sqrt(s2/n),n-1)./sqrt(s2/n);
plot(x,pm)
set(gca,'YTick',[],'XLim',[-45 45],'XTick',-40:20:40)
text(20,0.35,'\mu=26.2','FontSize',16,'FontWeight','bold')
xlim([-45 45])
title('Posterior of \mu')

pause

% currently accepted value
line([33 33],[0 0.3],'Color',[0 0.5 0])
text(35,0.33,['Currently     '; 'accepted value'],'HorizontalAlignment','left','FontSize',14,'FontWeight','bold')

pause

% plot again with extra plot 3
clf
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'Position',[600 60 700 640]);

% histogram of data
subplot(3,1,1)
set(gca,'box','on')
hist(y,-43:2:41)
set(gca,'YTick',[],'XTick',-40:20:40)
xlim([-45 45])
title('Newcomb''s measurements')

% the posterior distribution of mu
subplot(3,1,2)
set(gca,'box','on')
x=linspace(18.5,33.5,1000);
% the exact posterior distribution (t-distribution)
pm=tpdf((x-mean(y))/sqrt(s2/n),n-1)./sqrt(s2/n);
plot(x,pm)
set(gca,'YTick',[],'XLim',[-45 45],'XTick',-40:20:40)
title('Posterior of \mu')

% currently accepted value
line([33 33],[0 0.3],'Color',[0 0.5 0])
text(35,0.33,['Currently     '; 'accepted value'],'HorizontalAlignment','left','FontSize',14,'FontWeight','bold')

% What if we remove suspicious observations below zero?
% In real situations, use this only as approximation
% often there is better way of handling "outliers",
% for example, using more robust likelihood like t-distribution
y=y(y>0);
n=length(y);
s2=var(y);
my=mean(y);

% the posterior distribution of mu given y>0
subplot(3,1,3)
set(gca,'box','on')
x=linspace(18.5,33.5,1000);
% the exact posterior distribution (t-distribution)
pm=tpdf((x-mean(y))/sqrt(s2/n),n-1)./sqrt(s2/n);
plot(x,pm)
set(gca,'YTick',[],'XLim',[-45 45],'XTick',-40:20:40)
title('Posterior of \mu given y>0')

line([33 33],[0 0.63],'Color',[0 0.5 0])
text(35,0.33,['Currently     '; 'accepted value'],'HorizontalAlignment','left','FontSize',14,'FontWeight','bold')

% Note that even after removing "outliers" it seems that there
% was a systematic measurement error in Newcombs' experiment
