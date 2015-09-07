% Bayesian data analysis
% Aki Vehtari <Aki.Vehtari@aalto.fi>

% 437 girls and 543 boys have been observed. Calculate and plot the
% posterior distribution of the proportion of girls $\theta $, using
% uniform prior on $\theta $

% Set general figure properties
clf
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')
set(gcf,'DefaultLineLineWidth',2)

% Plot posterior Beta(438,544)
x=0.375:0.001:0.525;
a=438;b=544;
% betapdf computes the posterior density
pd=betapdf(x,a,b); 
% plot
plot(x,pd)
set(gca,'XLim',[0.375 0.525],'XTick',[0.4 0.45 0.485 0.5])
set(gca,'Ytick',[]);

% Plot the proportion of girl babies in general population
line([0.485 0.485],ylim,'LineStyle','--','Color','k')

% Colorize 95% posterior interval
hold on % hold previous plot
% the following line creates 100 evenly spaced values
% from 2.5% quantile to 97.5% quantile (i.e., 95% central interval)
% betainv computes the value for a given quantile given parameters a and b
xx=linspace(betainv(0.025,a,b),betainv(0.975,a,b),100);
% compute the posterior density
ppd=betapdf(xx,a,b);
% plot area
h=area(xx,ppd);
% add a legend
legend(h,'95% posterior interval')
% remove hold
hold off
% additional decoration
set(gca,'Ytick',[]);
title('Uniform prior -> Posterior is Beta(438,544)')

