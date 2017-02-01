% Bayesian data analysis
% Aki Vehtari <Aki.Vehtari@aalto.fi>

% Probability of a girl birth given placenta previa (BDA3 p. 37).
% Illustrate the effect of prior.
% Comparison of posterior distributions with different
% parameter values for the beta prior distribution. 
    
% Set general figure properties
clf
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')
set(gcf,'DefaultLineLineWidth',2)
x=0.375:0.001:0.525;

% Use 3 subplots
subplot(3,1,1)
% Plot posterior, data (437,543) and uniform prior (Beta(1,1))
a=438;b=544;
plot(x,betapdf(x,a,b))
% Plot prior and posterior: data (437,543) and prior Beta(0.485*2,(1-0.485)*2)
ap=0.485*2;bp=(1-0.485)*2;
aa=437+ap;bb=543+bp;
hold on
plot(x,betapdf(x,ap,bp),'k:',x,betapdf(x,aa,bb),'r')
hold off
h=legend('Post with unif prior','Informative prior','Posterior',2);
set(h,'FontSize',12,'FontWeight','normal')
% Additional decoration
set(gca,'XLim',[0.375 0.525],'XTick',[0.4 0.45 0.485 0.5])
set(gca,'Ytick',[]);
line([0.485 0.485],ylim,'LineStyle','--','Color','k')
title('alpha/(alpha+beta)=0.485, alpha+beta=2')
% wait for enter
pause

% Second subplot
subplot(3,1,2)
% Plot posterior, data (437,543) and uniform prior (Beta(1,1))
a=438;b=544;
plot(x,betapdf(x,a,b))
% Plot posterior, data (437,543) and prior Beta(0.485*20,(1-0.485)*20)
ap=0.485*20;bp=(1-0.485)*20;
aa=437+ap;bb=543+bp;
hold on
plot(x,betapdf(x,ap,bp),'k:',x,betapdf(x,aa,bb),'r')
hold off
% Additional decoration
set(gca,'XLim',[0.375 0.525],'XTick',[0.4 0.45 0.485 0.5])
set(gca,'Ytick',[]);
line([0.485 0.485],ylim,'LineStyle','--','Color','k')
title('alpha/(alpha+beta)=0.485, alpha+beta=20')
% wait
pause

% Third subplot
subplot(3,1,3)
% Plot posterior, data (437,543) and uniform prior (Beta(1,1))
a=438;b=544;
plot(x,betapdf(x,a,b))
% Plot posterior, data (437,543) and prior Beta(0.485*20,(1-0.485)*20)
ap=0.485*200;bp=(1-0.485)*200;
aa=437+ap;bb=543+bp;
hold on
plot(x,betapdf(x,ap,bp),'k:',x,betapdf(x,aa,bb),'r')
hold off
% Additional decoration
set(gca,'XLim',[0.375 0.525],'XTick',[0.4 0.45 0.485 0.5])
set(gca,'Ytick',[]);
line([0.485 0.485],ylim,'LineStyle','--','Color','k')
title('alpha/(alpha+beta)=0.485, alpha+beta=200')
