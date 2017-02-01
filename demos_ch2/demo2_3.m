% Bayesian data analysis
% Aki Vehtari <Aki.Vehtari@aalto.fi>

% Probability of a girl birth given placenta previa (BDA3 p. 37).
% Simulate samples from Beta(438,544), draw a histogram with
% quantiles, and do the same for a transformed variable.

% Set general figure properties
clf
set(gcf,'DefaultTextFontSize',16,'DefaultTextFontWeight','bold')
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')
set(gcf,'DefaultLineLineWidth',2)
x=0.375:0.001:0.525;

% Sample from posterior Beta(438,544)
% Get all samples at once and store them in vector 'r'
a=438;b=544;
th=betarnd(a,b,1,10000);

% Animate sampling from posterior
subplot(1,1,1)
% Plot posterior density
plot(x,betapdf(x,a,b))
% Additional decoration
set(gca,'XLim',[0.375 0.525],'XTick',[0.4 0.45 0.485 0.5])
set(gca,'Ytick',[]);
% Hold plot
hold on
% Plot 200 first samples
for i1=1:200
  h=line([th(i1) th(i1)],[0 betapdf(th(i1),a,b)],'LineStyle','-','Color','b');
  plot(th(i1),0,'b.')
  % Begin with a longer pause and then accelerate gradually
  pause(1.5/(1+i1))
  delete(h)
end
% Wait for enter
pause

% Plot the histogram of all samples of theta
subplot(2,1,1)
% hist produces histogram plot, second argument tells how many bins
hist(th,30)
% compute 2.5% and 97.5% quantile approximation using samples
th25=prctile(th,2.5);
th975=prctile(th,97.5);
% plot lines, ylim return limits of y axis
line([th25 th25],ylim)
line([th975 th975],ylim)
% decoration
xlabel('theta')
set(gca,'Ytick',[]);
yl=ylim;
text(th25,yl(2)+12,'2.5%','HorizontalAlignment','center')
text(th975,yl(2)+15,'97.5%','HorizontalAlignment','center')

% Plot histogram for phi=(1-theta)/theta (with all samples)
pause
subplot(2,1,2)
% hist produces histogram plot, second argument tells how many bins
phi=(1-th)./th;
hist(phi,30)
% compute 2.5% and 97.5% quantile approximation using samples
phi25=prctile(phi,2.5);
phi975=prctile(phi,97.5);
% plot lines, ylim return limits of y axis
line([phi25 phi25],ylim)
line([phi975 phi975],ylim)
% decoration
xlabel('phi=(1-theta)/theta')
set(gca,'Ytick',[]);
yl=ylim;
text(phi25,yl(2)+12,'2.5%','HorizontalAlignment','center')
text(phi975,yl(2)+15,'97.5%','HorizontalAlignment','center')
