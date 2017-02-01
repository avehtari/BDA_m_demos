% Bayesian data analysis
% Aki Vehtari <Aki.Vehtari@aalto.fi>

% Probability of a girl birth given placenta previa (BDA3 p. 37).
% Calculate the posterior distribution on a discrete grid of points by
% multiplying the likelihood and a non-conjugate prior at each point,
% and normalizing over the points. Simulate samples from the resulting
% non-standard posterior distribution using inverse cdf using the
% discrete grid.

% Set general figure properties
clf
set(gcf,'DefaultTextFontSize',16,'DefaultTextFontWeight','bold')
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')
set(gcf,'DefaultLineLineWidth',2)

% First subplot
subplot(3,1,1)
% Plot posterior with uniform prior
% data (437,543) and uniform prior (Beta(1,1))
a=438;b=544;
x=0.35:0.001:0.6;
plot(x,betapdf(x,a,b))
% decorations
xlim([0.35 0.6])
set(gca,'Ytick',[]);
text(0.36,max(ylim)*0.8,'Poster with uniform prior is Beta(438,544)')

% Pause for enter
pause
% Second subplot
subplot(3,1,2)
% Compute the density of non-conjugate prior in discrete points, i.e. in a grid
% this non-conjugate prior is the same as in figure 2.4 in the book
px=0:0.001:1;pp=ones(size(px));
pi1=find(px==0.385);
pi2=find(px==0.485);
pi3=find(px==0.585);
pm=11;
pp(pi1:pi2)=linspace(1,pm,length(pi1:pi2));
pp(pi3:-1:pi2)=linspace(1,pm,length(pi1:pi2));
% normalize the prior
pp=pp./sum(pp);
% plot non-conjugate prior (XAxis [0,1])
plot(px,pp)
% decorations
set(gca,'Ytick',[]);
text(0.03,max(ylim)*0.8,'Non-conjugate prior')

% pause for enter
pause
% plot non-conjugate prior again with tighter xlim (XAxis [0.35 0.6])
plot(px,pp)
% decorations
xlim([0.35 0.6])
set(gca,'Ytick',[]);
text(0.36,max(ylim)*0.8,'Non-conjugate prior')

% pause for enter
pause
% third subplot
subplot(3,1,3)
% compute the un-normalized non-conjugate posterior in a grid
po=betapdf(px,a,b).*pp;
% normalize the posterior
po=po./sum(po);
% plot the non-conjugate posterior (XAxis [0.35 0.6])
plot(px,po)
% decorations
xlim([0.35 0.6])
set(gca,'Ytick',[]);
text(0.36,max(ylim)*0.8,'Non-conjugate posterior')

pause

% next demonstrate inverse cdf sampling
% clear figure
clf
% first subplot
subplot(3,1,1)
% plot the non-conjugate posterior (XAxis [0.35 0.6])
plot(px,po)
% decorations
set(gca,'XLim',[0.35 0.6],'Xtick',[0.4:0.1:0.6])
set(gca,'Ytick',[]);
text(0.47,max(ylim)*0.8,'Non-conjugate posterior')
title('Illustration of inverse-cdf-method')

pause

% second subplot
subplot(3,1,2)
% compute the cumulative density in a grid
pc=cumsum(po);
% Plot non-conjugate posterior-cdf (XAxis [0.35 0.6])
plot(px,pc)
% decorations
set(gca,'XLim',[0.35 0.6],'Xtick',[0.4:0.1:0.6])
set(gca,'Ytick',[0 0.5 1]);
h=text(0.47,max(ylim)*0.8,'Posterior-cdf');
% draw now
drawnow

% Inverse-cdf sampling
% In matlab it is more efficient to call rand once and get vector
% of random numbers
% rand returns uniform random numbers [0,1]
% rng(state) is used to set state of the randon number generator
rng(2601)
r=rand(1,10000);
% for inverse-cdf we need to start cdf from exact zero
% for plotting this did not matter
qc=[0 pc];
clear rr
for i1=1:10000
  % finding the value for which cdf is equal to the sampled uniform number
  rr(i1)=px(sum(qc<r(i1)));
end
% alternatively use catrand which is usually faster approach
% rr=px(catrand(po,1,10000));

% wait for enter
pause

% now animate inverse-cdf sampling
% start with illustrating first iteration
i1=1;
h1=line([0 rr(i1)],[r(i1) r(i1)],'LineStyle','--','Color','b');
ht=text(0.352,r(i1),'Sample uniformly from [0,1]','FontSize',14,'Fontweight','bold','VerticalAlign','top');
% pause
pause
delete(ht)
h2=line([rr(i1) rr(i1)],[r(i1) 0],'LineStyle','--','Color','b');
line(rr(i1),0,'LineStyle','none','Color','b','Marker','.')
ht=text(rr(i1),0.1,' Find corresponding value from cdf','FontSize',14,'Fontweight','bold');
% wait for enter
pause

% delete the first illustration
delete(h1)
delete(h2)
delete(ht)

% animate inverse cdf-sampling for next 99 iterations
for i1=2:100
  h1=line([0 rr(i1)],[r(i1) r(i1)],'LineStyle','--','Color','b');
  h2=line([rr(i1) rr(i1)],[r(i1) 0],'LineStyle','--','Color','b');
  line(rr(i1),0,'LineStyle','none','Color','b','Marker','.')
  % pause for a decreasing amount of time
  pause(1.5/(1+i1))
  delete(h1)
  delete(h2)
end

pause

% third subplot
subplot(3,1,3)
% and finally plot histogram of posterior samples
% data (437,543) and non-conjugate prior
hist(rr,30)
% decorations
xlim([0.35 0.6])
set(gca,'Ytick',[]);
h=text(0.47,max(ylim)*0.8,'Histogram of posterior samples');
