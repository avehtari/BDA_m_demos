% Bayesian data analysis
% Aki Vehtari <Aki.Vehtari@aalto.fi>

% Hierarchical model for SAT-example data (BDA3, p. 102)

% SAT-example data (BDA3 p. 120)
% y is the estimated treatment effect
% s is the standard error of effect estimate
y=[28  8  -3  7  -1  1  18  12];
s=[15  10  16  11  9  11  10  18];

clf
% Default settings for figure
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'Position',[200 60 820 640]);
box on

% the number of observations
M=length(y);

% Plot data, use normpdf to plot both the y_j and sigma_j
clf
x=linspace(-60,80,1000);
box on
for i=1:M
  px=normpdf(x,y(i),s(i));
  line(x,px)
end
set(gca,'Ytick',[])
title('Data')
xlabel('Treatment effect')
line([0 0],ylim,'LineStyle',':','Color','k','LineWidth',.5)

% Plot separate, pooled, and hierarchical model
pause
clf
subplot(3,1,1)
x=linspace(-40,60,1000);
box on
for i=1:M
  px=normpdf(x,y(i),s(i));
  h1=line(x,px);
end
i=1;
px=normpdf(x,y(i),s(i));
h2=line(x,px,'color','r');
set(gca,'Ytick',[])
xlim([-40 60])
title('Separate model')
hl=legend([h2 h1],'School A','Other schools');
set(hl,'FontSize',14)


subplot(3,1,2)
box on
px=normpdf(x,sum(y./s.^2)/sum(1./s.^2),sqrt(1./sum(1./s.^2)));
line(x,px)
h2=line(x,px,'color','r');
set(gca,'Ytick',[])
title('Pooled model')
xlim([-40 60])
hl=legend(h2,'All schools');
set(hl,'FontSize',14)

subplot(3,1,3)
x=linspace(-40,60,500);
box on
% Load the pre-computed results for the hierarchical model
% Replace this with your own code in the related exercise
load demo5_2
%>> whos -file demo5_2
%  Name      Size                    Bytes  Class
%  pxm       8x500                   32000  double array
%  t         1x1000                   8000  double array
%  tm        8x1000                  64000  double array
%  tp        1x1000                   8000  double array
%  tsd       8x1000                  64000  double array
for i=1:M
  h1=line(x,pxm(i,:));
end
h2=line(x,pxm(1,:),'Color','r');
set(gca,'Ytick',[])
title('Hierarchical model')
hl=legend([h2 h1],'School A','Other schools');
set(hl,'FontSize',14)

% Plot various marginal and conditional posterior summaries
pause
clf
subplot(3,1,1)
plot(t,tp)
set(gca,'Ytick',[])
title('Marginal posterior density p(\tau|y)')
xlim([0 35])
xlabel('\tau')
ylabel('p(\tau|y)')

subplot(3,1,2)
h1=plot(t,tm,'b',t,tm(1,:),'r');
xlim([0 35])
title('Conditional posterior means of effects E(\theta_j|\tau,y)')
hl=legend([h1(end) h1(1)],'School A','Other schools','location','northwest');
set(hl,'FontSize',14)
xlabel('\tau')
ylabel('E(\theta_j|\tau,y)')

subplot(3,1,3)
h1=plot(t,tsd,'b',t,tsd(1,:),'r');
xlim([0 35])
title('Conditional posterior standard deviations of effects sd(\theta_j|\tau,y)')
hl=legend([h1(end) h1(1)],'School A','Other schools','location','northwest');
set(hl,'FontSize',14)
xlabel('\tau')
ylabel('std(\theta_j|\tau,y)')
