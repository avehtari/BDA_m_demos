% Bayesian data analysis
% Aki Vehtari <Aki.Vehtari@aalto.fi>

% PSRF demonstration

clf
% default settings for figure
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')

y1=0;y2=0;r=0.8; % parameters of the toy posterior (normal-distribution), r=rho
S=[1 r;r 1];     % covariance matrix the toy posterior (normal-distribution)
t1=-2.5;t2=2.5;  % starting value for the chain
sp=0.8;          % scale (std) for proposal distribution (use normal distribution)
M=10000;         % number of iterations
tt=zeros(M,2);   % variable for saving samples, tt contains joint samples for theta_1 and theta_2
tt(1,:)=[t1 t2]; % save starting point
rr=0;            % variable for saving number of rejections


% 10 different starting points (t1=theta_1, t2=theta_2)
t1s=[-2.5 -2.5 -2.5 0   -1  1  0   2.5 2.5 2.5];
t2s=[ 2.5  0   -2.5 2.5  1 -1 -2.5 2.5 0  -2.5];

%
clf
% Plot 90% HPD
% see BDA3 p. 85, for how to compute HPD for multivariate normal
% in 2d-case contour for 90% HPD is an ellipse, whose semimajor
% axes can be computed from the eigenvalues of the covariance
% matrix scaled by a value selected to get ellipse match the
% density at the edge of 90% HPD. Angle of the ellipse could be 
% computed from the eigenvectors, but since the marginals are same
% we know that angle is pi/4
q=sort(sqrt(eig(S))*2.147);
he=ellipse(q(2),q(1),pi/4,0,0,'r'); 
% Set axes
set(gca,'DataAspectRatio',[1 1 1])
set(gca,'XLim',[-4.5 4.5],'YLim',[-4.5 4.5])
set(gca,'Box','on')
% plot the starting points
h1=line(t1s,t2s,'Marker','o','LineStyle','none');
title('Convergence demonstration')
legend([he h1],'90% HPD','Starting points')

% load pre-run Metropolis chains
% note that proposal distribution was intentionally selected to be
% slightly too small, to better illustrate convergence diagonstics
% Since, implementation of the Metropolis algorithm is on of the
% exercises, we load here pre-computed chains and PSRF-values

% `tts' contains samples, `p1' and 'p2' contain PSRF values for t1
% and t2 using 50% burn-in, `pp1' and 'pp2' contain PSRF values for
% t1 and t2 using 10% burn-in, PSRF-values have been computed for
% each time-step. 
load demo11_4
%>> whos -file demo11_4
%  Name      Size                           Bytes  Class
%  p1        1x10000                        80000  double array
%  p2        1x10000                        80000  double array
%  pp1       1x10000                        80000  double array
%  pp2       1x10000                        80000  double array
%  tts   20000x2x10                       3200000  double array
% PSRF values have been pre-computed because computing them for
% each time step takes some time. Computation used is shown below
% for n2=10:10:M
%  n1=n2/2;
%  p1(n2)=psrf(tts(n1:n2,1,:));
%  p2(n2)=psrf(tts(n1:n2,2,:));
% end
% for n2=10:10:M
%  n1=n2/10;
%  pp1(n2)=psrf(tts(n1:n2,1,:));
%  pp2(n2)=psrf(tts(n1:n2,2,:));
% end
pause

% plot 50 first iterations from the chains
for ti=2:50
  for i1=1:10
    line(tts(ti-1:ti,1,i1),tts(ti-1:ti,2,i1))
  end
  pause(1/ti)
end
% no convergence yet
title('Convergence demonstration - not yet converged')
pause

% Plot trends of 50 first samples
subplot(2,1,1)
plot(squeeze(tts(1:50,1,:)))
ylabel('\theta_1')
title('Convergence demonstration - not yet converged')
subplot(2,1,2)
plot(squeeze(tts(1:50,2,:)))
ylabel('\theta_2')
% no convergence yet
pause

% Plot trends from 500th sample to end
set(gcf,'DefaultLineLineWidth',1)
subplot(2,1,1)
plot(501:M,squeeze(tts(501:M,1,:)))
ylabel('\theta_1')
title('Convergence demonstration - visually converged')
subplot(2,1,2)
plot(501:M,squeeze(tts(501:M,2,:)))
ylabel('\theta_2')
set(gcf,'DefaultLineLineWidth',2)
% seems to be converged
pause

% Plot PSRF for different chain lengths
% remove first half as a burn-in
subplot(2,1,1)
semilogx(10:10:M,p1(10:10:M))  
axis([10 M 0 30])
title('Running PSRF')
ylabel('\theta_1')
legend('PSRF(n/2:n)')
subplot(2,1,2)
semilogx(10:10:M,p2(10:10:M))  
axis([10 M 0 30])
ylabel('\theta_2')
legend('PSRF(n/2:n)')

pause
% change of axis
subplot(2,1,1)
hl=line([10 M],[1 1],'LineStyle','--','Color','k');
axis([10 M 0 5])
subplot(2,1,2)
hl=line([10 M],[1 1],'LineStyle','--','Color','k');
axis([10 M 0 5])

pause
% change of axis
subplot(2,1,1)
axis([1000 M .5 1.5])
subplot(2,1,2)
axis([1000 M .5 1.5])
pause

% plot PSRF with 50% burn-in and 10% burn-in
subplot(2,1,1)
plot(10:10:M,p1(10:10:M),10:10:M,pp1(10:10:M))  
hl=line([4 M],[1 1],'LineStyle','--','Color','k');
axis([200 M 0.9 1.25])
title('Running PSRF with different burn-in length')
ylabel('\theta_1')
legend('PSRF(n/2:n)','PSRF(n/10:n)')
subplot(2,1,2)
plot(10:10:M,p1(10:10:M),10:10:M,pp2(10:10:M))  
hl=line([4 M],[1 1],'LineStyle','--','Color','k');
axis([200 M 0.9 1.25])
ylabel('\theta_2')
legend('PSRF(n/2:n)','PSRF(n/10:n)')
pause

% Illustrate the within and between variations
% and their progress during sampling
clf
Y1=linspace(-4,4,200);
subplot(2,1,1)
for n2=[10 20 50 100 200 500 1000 2000 5000 10000]
  clf
  n1=n2/2;
  [R,neff,V,W,B]=psrf(tts(n1:n2,1,:));
  qt=squeeze(tts(n1:n2,1,:));
  subplot(2,1,1)
  title('Illustration of within and between variations')
  for i1=1:10
    h1=line(Y1,normpdf(Y1,mean(qt(:,i1)),std(qt(:,i1))));
  end
  legend(h1,'Single chains')
  text(-3.8,max(ylim)*0.9,sprintf('m*n=%d, m*n/2=%d, neff=%d',n2*10,n2*10/2,floor(neff)),...
       'FontSize',16,'FontWeight','bold')
  qm=mean(qt(:));
  subplot(2,1,2)
  h2=line(Y1,normpdf(Y1,qm,sqrt(V)),'Color',[0 0.5 0]);
  h3=line(Y1,normpdf(Y1,qm,sqrt(W)),'Color','r');
  legend([h2 h3],'V','W')
  text(-3.8,max(ylim)*0.9,sprintf('PSRF=%.2f',R),'FontSize',16,'FontWeight','bold')
  pause
end

