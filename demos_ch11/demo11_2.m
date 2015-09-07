% Bayesian data analysis
% Aki Vehtari <Aki.Vehtari@aalto.fi>

% Metropolis sampling demonstration

% new figure
figure
% Default settings for figure
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')

y1=0;y2=0;r=0.8; % parameters of the Normal-distribution, which is used as a
                 % toy target distribution for illustration
S=[1 r;r 1];     % covariance matrix of the toy distribution

t1=-2.5;t2=2.5;  % starting value for the chain
sp=0.8;          % scale for the proposal distribution
M=5000;          % number of iterations
                 % unlike in my Gibbs implementation, here all parameters
                 % are updated at once, and thus one iteration is complete
                 % iteration
tt=zeros(M,2);   % variable for saving samples
tt(1,:)=[t1 t2]; % save starting point
rr=0;            % variable for saving number of rejections

% Illustrate the target distribution
clf
Y1=linspace(-4.5,4.5,200);
Y2=linspace(-4.5,4.5,200);
lp=zeros(200);
for i1=1:length(Y1)
  for i2=1:length(Y2)
    lp(i1,i2)=mnorm_lpdf([Y1(i1) Y2(i2)],[y1 y2],S);
  end
end
pcolor(Y1,Y2,exp(lp))
shading interp
set(gca,'DataAspectRatio',[1 1 1])
set(gca,'Box','on')
set(gcf,'renderer','painters')
xlabel('\theta_1')
ylabel('\theta_2')
title('Toy target distribution')
pause

% Metropolis sampling here
% For demonstration load pre-computed values
% Replace this with your algorithm!
% `ttall' is all proposals, `all' is all acceptance ratios
load demo11_1
%whos -file demo11_1
%  Name        Size                    Bytes  Class
%  all         1x4999                  39992  double array
%  tt       5000x2                     80000  double array
%  ttall    4999x2                     79984  double array

% rejection rate
% your algorithm should update number of rejections 'rr'
% divide the number of rejections with M to get rejection rate
%rr/M

% The rest is just for illusttration

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

% plot the starting point
hs=line(tt(1,1),tt(1,2),'Marker','o','LineStyle','none');
title('Metropolis-algorithm')
xlabel('\theta_1')
ylabel('\theta_2')
hl=legend([he hs],'90% HPD','Starting point',3);
set(hl,'FontSize',14)
pause

% plot progession during 20 first iterations
samplestring='Starting point';
ht=text(-4,4,'','FontSize',16,'FontWeight','bold','Color','b');
for i1=1:20
  % plot the 90% HPD of the proposal distribution
  h1=ellipse(sp*2.147,sp*2.147,0,tt(i1,1),tt(i1,2),'b');
  hl=legend([he hs h1],'90% HPD',samplestring,'Proposal distribution (90% HPD)',3);
  set(hl,'FontSize',14)
  pause
  % plot the proposal and ratio
  h2=line(ttall(i1,1),ttall(i1,2),'LineStyle','none','Marker','*');
  set(ht,'String',sprintf('Step %d - Acceptance ratio r=%.2f',i1,all(i1)),'Color','b')
  hl=legend([he hs h1 h2],'90% HPD',samplestring,'Proposal distribution (90% HPD)','Proposal \theta^*',3);
  set(hl,'FontSize',14)
  pause
  if ttall(i1,1)==tt(i1+1,1) && ttall(i1,2)==tt(i1+1,2)
    % accepted
    set(h2,'Color',[0 0.5 0]);
    hl=legend([he hs h1 h2],'90% HPD',samplestring,'Proposal distribution (90% HPD)','Proposal \theta^* accepted',3);
    set(hl,'FontSize',14)
    set(ht,'String',sprintf('Step %d - Proposal accepted, r=%.2f',i1,all(i1)),'Color',[0 0.5 0])
  else
    % rejected
    set(h2,'Color','r');
    hl=legend([he hs h1 h2],'90% HPD',samplestring,'Proposal distribution (90% HPD)','Proposal \theta^* rejected',3);
    set(hl,'FontSize',14)
    set(ht,'String',sprintf('Step %d - Proposal rejected, r=%.2f',i1,all(i1)),'Color','r')
  end
  pause
  delete(h1)
  delete(h2)
  h2=line(tt(i1+1,1),tt(i1+1,2),'LineStyle','none','Marker','o','Color','b');
  samplestring='Samples from the chain';
  hl=legend([he h2],'90% HPD',samplestring,3);
  set(hl,'FontSize',14)
  set(ht,'String',sprintf('Step %d',i1),'Color','b')
  pause
end

pause

% connect the samples with a line
h1=line(tt(1,1),tt(1,2),'Marker','o');
line(tt(1:21,1),tt(1:21,2));
hl=legend([he h1],'90% HPD','Markov chain',3);
set(hl,'FontSize',14)

pause

% plot progession up to 200 first iterations
for i1=22:201
  line(tt(i1-1:i1,1),tt(i1-1:i1,2));
  h2=line(tt(i1,1),tt(i1,2),'Marker','o','LineStyle','none');
  if ttall(i1,1)==tt(i1+1,1) && ttall(i1,2)==tt(i1+1,2)
    % accepted
    set(h2,'Color',[0 0.5 0]);
    set(ht,'String',sprintf('Step %d',i1),'Color',[0 0.5 0])
  else
    % rejected
    set(h2,'Color','r');
    set(ht,'String',sprintf('Step %d',i1),'Color','r')
  end
  set(h2,'Color','b');
  pause(2/(0.5+i1));
end
% Replace 53% with your rejection rate
set(ht,'Color','b');
text(-1.5,4,'Rejection rate 53%','FontSize',16,'FontWeight','bold');

pause
% Clear
clf
% Plot 90% HPD
q=sort(sqrt(eig(S))*2.147);
he=ellipse(q(2),q(1),pi/4,0,0,'r'); 
% Set axes
set(gca,'DataAspectRatio',[1 1 1])
set(gca,'XLim',[-4.5 4.5],'YLim',[-4.5 4.5])
set(gca,'Box','on')
% Plot 200 first samples without lines
h=line(tt(1:201,1),tt(1:201,2),'Marker','o','LineStyle','none');
title('Metropolis-algorithm')
hl=legend([he h],'90% HPD','Samples from the chain',3);
set(hl,'FontSize',14)

pause
% remove plot of the 200 first samples
delete(h)
% Remove burn-in and plot the every 5th sample starting from iteration 105
h=line(tt(105:5:end,1),tt(105:5:end,2),'Marker','o','LineStyle','none');
title('Metropolis-algorithm')
hl=legend([he h],'90% HPD','Samples from the chain after 200 step burn-in',3);
set(hl,'FontSize',14)

pause
% Clear
clf

subplot(3,1,1)
% show trends
burnin=105;
h1=plot([burnin:M],tt(burnin:M,:));
line([burnin M],[0 0])
axis tight
title('Metropolis-algorithm - trends')
hl=legend(h1,'\theta_1', '\theta_2');
set(hl,'FontSize',14)

subplot(3,1,2)
% show how the average behaves when more samples are obtained
burnin=105;
h1=plot([burnin:M],cumsum(tt(burnin:M,:))./repmat([1:length(burnin:M)]',1,2));
line([burnin M],[0 0])
axis([105 5000 -0.5 0.2])
title('Cumulative average')
hl=legend(h1,'\theta_1', '\theta_2',4);
set(hl,'FontSize',14)

subplot(3,1,3)
% show the autocorrelation
maxlag=20;
ac=acorr(tt(burnin:2:M,:),maxlag);
% Note: autocorrelation for 0 is 1
h1=plot(0:maxlag,[1 1; ac]);
set(gca,'XLim',[0 maxlag],'YLim',[-0.2 1])
line([0 maxlag],[0 0],'LineStyle','--','Color','k')
title('Estimate of the autocorrelation function')
hl=legend(h1,'\theta_1', '\theta_2');
set(hl,'FontSize',14)

pause
clf
% show how the average behaves when more samples are obtained
burnin=105;
h=plot([burnin:M],cumsum(tt(burnin:M,:))./repmat([1:length(burnin:M)]',1,2));
line([burnin M],[0 0])
axis([105 5000 -0.5 0.5])
legend('\theta_1','\theta_2')
title('Metropolis sampling - cumulative average')
% compute the effective number of samples
[R,neff] = psrf(tt(burnin:M,:));
% for a simpler plot, use the minimum neff (each parameter has its own neff)
neff=min(neff);
% plot approximative 95% quantiles for Markov chain Monte Carlo error
% this takes into account that samples are not independent
line([burnin:M],1.96./sqrt([burnin:M]/(M/neff)),'LineStyle','--','Color','k')
h2=line([burnin:M],-1.96./sqrt([burnin:M]/(M/neff)),'LineStyle','--','Color','k');
legend([h' h2],'\theta_1','\theta_2','95% interval for MCMC error')
pause
% For comparison plot 95% quantiles for Monte Carlo error
% if we would be able to obtain independent samples
line([burnin:M],1.96./sqrt([burnin:M]),'LineStyle','-.','Color','k')
h3=line([burnin:M],-1.96./sqrt([burnin:M]),'LineStyle','-.','Color','k');
legend([h' h2 h3],'\theta_1','\theta_2','95% interval for MCMC error','95% interval for independent MC')
