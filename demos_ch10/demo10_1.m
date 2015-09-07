% Bayesian data analysis
% Aki Vehtari <Aki.Vehtari@aalto.fi>

% Rejection sampling example

% a fake interesting distribution
% named q, to emphesize that it does not need to be normalized
x=linspace(-3,3,200);
r=[1.1 1.3 -0.1 -0.7 0.2 -0.4 0.06 -1.7 1.7 0.3 0.7 1.6 -2.06 -0.74 0.2 0.5]';
q=kernelp(r,x);

% Rejection sampling example
clf
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')
x=linspace(-3,3,200);
g=normpdf(x);
% M is computed by discrete approximation
M=max(q./g);
% prescale
g=M*g;
plot(x,q,x,g,'--')
hold on;
area(x,q,'FaceColor',[0.7 0.7 1]);
hold off
set(gca,'YTick',[])
legend('{ q(\theta|y)}','{ Mg(\theta)}')
title('Rejection sampling')
r1=-0.8;
[z,zi]=min(abs(x-r1));
h1=line([r1 r1],[0 q(zi)]);
h2=line([r1 r1],[q(zi) g(zi)],'LineStyle','--');
h3=text(r1,q(zi),'\leftarrow { q(\theta=r|y)}','FontSize',16,'FontWeight','bold');
h4=text(r1,g(zi),'\leftarrow { Mg(\theta=r)}','FontSize',16,'FontWeight','bold');
r2=0.3*g(zi);
line(r1,r2,'Marker','o','MarkerSize',8,'MarkerFaceColor',[0 0.5 0],'MarkerEdgeColor',[0 0.5 0])

pause
plot(x,q,x,g,'--')
hold on
area(x,q,'FaceColor',[0.7 0.7 1]);
hold off
set(gca,'YTick',[])
legend('{ q(\theta|y)}','{ Mg(\theta)}')
title('Rejection sampling')
r1=-0.8;
[z,zi]=min(abs(x-r1));
h1=line([r1 r1],[0 q(zi)]);
h2=line([r1 r1],[q(zi) g(zi)],'LineStyle','--');
h3=text(r1,q(zi),'\leftarrow { q(\theta=r|y)}','FontSize',16,'FontWeight','bold');
h4=text(r1,g(zi),'\leftarrow { Mg(\theta=r)}','FontSize',16,'FontWeight','bold');
r2=0.8*g(zi);
line(r1,r2,'Marker','o','MarkerSize',8,'MarkerFaceColor','r','MarkerEdgeColor','r')

pause
plot(x,q,x,g,'--')
hold on
area(x,q,'FaceColor',[0.7 0.7 1]);
hold off
set(gca,'YTick',[],'XLim',[-3 3])
legend('{ q(\theta|y)}','{ Mg(\theta)}')
title('Rejection sampling')
for i1=1:200;
  r1(i1)=randn;
  [z,zi(i1)]=min(abs(x-r1(i1)));
  h1=line([r1(i1) r1(i1)],[0 q(zi(i1))]);
  h2=line([r1(i1) r1(i1)],[q(zi(i1)) g(zi(i1))],'LineStyle','--');
  r2(i1)=rand*g(zi(i1));
  if r2(i1)>q(zi(i1))
    h3(i1)=line(r1(i1),r2(i1),'Marker','o','MarkerSize',8,'MarkerFaceColor','r','MarkerEdgeColor','r');
  else
    h3(i1)=line(r1(i1),r2(i1),'Marker','o','MarkerSize',8,'MarkerFaceColor',[0 0.5 0],'MarkerEdgeColor',[0 0.5 0]);
  end
  pause(2/(0.5+i1))
  delete(h1)
  delete(h2)
end
pause
delete(h3(r2>q(zi)))
pause
h4=line(repmat(r1(r2<=q(zi)),2,1),[zeros(size(r1(r2<=q(zi)))); r2(r2<=q(zi))],'Color','b');
line(r1(r2<=q(zi)),zeros(size(r1(r2<=q(zi)))),'LineStyle','none','Marker','o','MarkerSize',6,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1])
pause
delete(h3(r2<=q(zi)))
delete(h4)

pause
clf
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')
x=linspace(-3,3,200);
g=ones(size(x))*max(q);
g(x<=-1.5)=linspace(q(1),max(q(x<=-1.5)),length(x(x<=-1.5)));
g(x>-1.5&x<=-0.1)=linspace(max(q(x<=-1.5)),max(q(x>-1.5&x<=-0.1)),length(x(x>-1.5&x<=-0.1)));
g(x>-0.1&x<=2.3)=linspace(max(q(x>-1.5&x<=-0.1)),max(q(x>2.3)),length(x(x>-0.1&x<=2.3)));
g(x>2.3)=linspace(max(q(x>2.3)),q(end),length(x(x>2.3)));
M=max(q./g);
plot(x,q,x,M*g,'--')
hold on
area(x,q,'FaceColor',[0.7 0.7 1]);
hold off
set(gca,'YTick',[],'XTick',[])
legend('q(\theta|y)','Mg(\theta)')
title('Rejection sampling - alternative proposal distribution')

