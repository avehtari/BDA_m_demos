% Bayesian data analysis
% Aki Vehtari <Aki.Vehtari@aalto.fi>

% Posterior predictive checking
% Binomial example - Testing sequential dependence example

clf
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')

% Testing sequential dependence example (Gelman et al p. 163)
y=[1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0];
Ty=sum(diff(y)~=0);

% sufficient statistics
n=length(y);
s=sum(y);

fprintf('Wait')
for i=1:10000
  t=betarnd(s+1,n-s+1);
  yr=binornd(1,t,n,1);
  Tyr(i)=sum(diff(yr)~=0);
  if ~rem(i,1000)
    fprintf('.')
  end
end

clf
hist(Tyr,0:17)
xlim([-.5 17.5])
set(gca,'YTick',[])
h1=get(gca,'Children');
h2=line([Ty Ty],ylim);
legend([h2 h1],'T(y)','T(yrep)')
mean(Tyr<=Ty)

title({'Binomial example - number of changes?','Pr(T(yrep,\theta)\leq T(y,\theta)|y)=0.03'})
