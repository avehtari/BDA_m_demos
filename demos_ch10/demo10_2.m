% Bayesian data analysis
% Aki Vehtari <Aki.Vehtari@aalto.fi>

% Importance sampling example

% a fake interesting distribution
% named q, to emphesize that it does not need to be normalized
r=[1.1 1.3 -0.1 -0.7 0.2 -0.4 0.06 -1.7 1.7 0.3 0.7 1.6 -2.06 -0.74 0.2 0.5]';
x=linspace(-3,3,200);
q=kernelp(r,x);

% Importance sampling example
clf
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')
x=linspace(-3,3,200);
g=normpdf(x);
subplot(2,1,1)
plot(x,q,x,g,'--')
set(gca,'YTick',[])
legend('{q(\theta|y)}','{g(\theta)}')
title('Importance sampling')
subplot(2,1,2)
w=q./g;
plot(x,w,'-.')
legend('q(\theta|y)/g(\theta)')
title('Importance weight function')

pause
r=randn(1,100);
ri=binsgeq(x,r);
h=plot(x,w,'-.');
hh=line([x(ri); x(ri)],[zeros(1,100); w(ri)],'LineWidth',2,'Color','b');
xlim(xlim)
title('Importance weights')
pause
delete(h) 