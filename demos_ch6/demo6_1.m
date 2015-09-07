% Bayesian data analysis
% Aki Vehtari <Aki.Vehtari@aalto.fi>

% Posterior predictive checking demo
clf
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')

% Light speed example
y=load('light.txt');
% sufficient statistics
n=length(y);
s2=var(y);
s=sqrt(s2);
my=mean(y);

clf
% Create 10 random replicate data sets from the posterior predictive
% density
% Each set has same number of virtual observations as the original
% data set
for i=1:10
  % n samples from the posterior predicitve density
  pp=trnd(n-1,n,1).*sqrt(1+1/n).*s+my;
  subplot(5,2,i)
  set(gca,'box','on')
  hist(pp,-45:5:55)
  set(gca,'YTick',[])
  set(gca,'XLim',[-48 58])
end
% Replace one of the replicates with observed data
% If you can spot which one has been replaced, it means that the
% replicates do not resemble the original data and thus the model has
% a defect
i=ceil(rand*10);
subplot(5,2,i)
set(gca,'box','on')
hist(y,-45:5:55)
set(gca,'YTick',[])
set(gca,'XLim',[-48 58])
subplot(5,2,1)
h=title({'Light speed example: Observed data + Replicated datasets',...
         'Can you spot which one is the observed data?'});
set(h,'Position',get(h,'Position')+[45 0 0])
drawnow

% Generate 1000 replicate data sets
for i=1:1000
  pps(:,i)=trnd(n-1,n,1).*sqrt(1+1/n).*s+my;
end
fprintf('.')
pause
% Minimum value from each replicate data sets
pp=min(pps);
clf
% Plot the distribution of the minimum values (20 = number of bars in histogram)
hist(pp,20);
% Plot the minimum observed value in the original data set
h1=get(gca,'Children');
h2=line([min(y) min(y)],ylim);
legend([h2 h1],'The smallest observation of the original data','The smallest observations of the replicated datasets')

