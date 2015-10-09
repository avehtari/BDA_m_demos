% Bayesian data analysis
% Aki Vehtari <Aki.Vehtari@aalto.fi>

% Gibbs sampling demonstration

% new figure
figure
% Default settings for figure
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')

y1=0;y2=0;r=0.8; % parameters of the Normal-distribution, which is used as a
                 % toy target distribution for illustration
S=[1 r;r 1];     % covariance matrix of the toy distribution

t1=-2.5;t2=2.5;  % starting value for the chain
M=2*1000;        % number of iterations
                 % note that in my implementation one iteration updates
                 % only one parameter, and one complete iteration
                 % updating both parameters takes two basic iterations
                 % this implementation was used to make plotting of Gibbs
                 % sampler's zig-zagging
                 % in plots I have showed iteration/2 to show only
                 % complete iterations
                 % you can implement this also by saving only the final
                 % state of complete iteration updating all parameters
tt=zeros(M,2);   % variable for saving samples
tt(1,:)=[t1 t2]; % save starting point

tt=zeros(M,2);   % variable for saving samples
tt(1,:)=[t1 t2]; % save starting point

% these are used for plotting
Y1=linspace(-4.5,4.5,200);
Y2=linspace(-4.5,4.5,200);

% Gibbs sampling here
% For demonstration load pre-computed values
% Replace this with your algorithm!
% tt is a M x 2 matrix, with M samples of both theta_1 and theta_2
load demo11_2
%>> whos -file demo11_2
%  Name      Size                    Bytes  Class
%  tt     2001x2                     32016  double array

% The rest is just for illustration

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
h1=line(tt(1,1),tt(1,2),'Marker','o','LineStyle','none');
title('Gibbs sampling')
xlabel('\theta_1')
ylabel('\theta_2')
hl=legend([he h1],'90% HPD','Starting point','location','southeast');
set(hl,'FontSize',14)
pause
for i1=1:2:5
  % plot the conditional distribution of theta_1 given theta_2
  hl1=line(xlim,[tt(i1,2) tt(i1,2)]','Color','k','LineStyle','--','LineWidth',2);
  hl2=line(Y1,tt(i1,2)+normpdf(Y1,y1+r.*(tt(i1,2)-y2),sqrt((1-r.^2)))*1,'Color','b','LineStyle','-','LineWidth',2);
  if i1==1
    hl=legend([he h1 hl2],'90% HPD','Starting point','Conditional density given \theta_2','location','southeast');
  else
    hl=legend([he h1 hl2],'90% HPD','Samples from the chain','Conditional density given \theta_2','location','southeast');
  end
  set(hl,'FontSize',14)
  pause
  % sample
  line(tt(i1+1,1),tt(i1+1,2),'LineStyle','none','Marker','o','Color','b')
  hl=legend([he h1 hl2],'90% HPD','Samples from the chain','Conditional density given \theta_2','location','southeast');
  set(hl,'FontSize',14)
  pause
  delete(hl1)
  delete(hl2)
  % plot the conditional distribution of theta_2 given theta_1
  hl1=line([tt(i1+1,1) tt(i1+1,1)]',get(gca,'YLim'),'Color','k','LineStyle','--','LineWidth',2);
  hl2=line(tt(i1+1,1)+normpdf(Y2,y2+r.*(tt(i1+1,1)-y1),sqrt((1-r.^2)))*1,Y2,'Color','b','LineStyle','-','LineWidth',2);
  hl=legend([he h1 hl2],'90% HPD','Samples from the chain','Conditional density given \theta_1','location','southeast');
  set(hl,'FontSize',14)
  pause
  % plot the sample
  line(tt(i1+2,1),tt(i1+2,2),'LineStyle','none','Marker','o','Color','b')
  pause
  delete(hl1)
  delete(hl2)
end
hl=legend([he h1],'90% HPD','Samples from the chain','location','southeast');
set(hl,'FontSize',14)

pause

% connect plotted samples with a line
h1=line(tt(1,1),tt(1,2),'Marker','o');
line(tt(1:7,1),tt(1:7,2));
hl=legend([he h1],'90% HPD','Markov chain','location','southeast');
set(hl,'FontSize',14)

pause

% Clear
clf
% Plot 90% HPD 
q=sort(sqrt(eig(S))*2.147);
he=ellipse(q(2),q(1),pi/4,0,0,'r'); 
% set axes
set(gca,'DataAspectRatio',[1 1 1])
set(gca,'XLim',[-4.5 4.5],'YLim',[-4.5 4.5])
set(gca,'Box','on')
% Plot again the first 7 samples
h1=line(tt(1,1),tt(1,2),'Marker','o'); % this is to get nice legend
line(tt(1:2:7,1),tt(1:2:7,2),'Marker','o','LineStyle','none');
line(tt(1:7,1),tt(1:7,2));
title('Gibbs sampling')
hl=legend([he h1],'90% HPD','Markov chain','location','southeast');
set(hl,'FontSize',14)

% Plot the next about 200 samples
% note that we could count the samples after each single update or
% after updating all the parameters, here we have stored for the illustration
% the samples after each single update but illustrate as we would have
% stored after updating all the parameters
ht=text(-4,4,sprintf('Step %d',(7-1)/2),'FontSize',16,'FontWeight','bold');
for i1=8:201
  line(tt(i1-1:i1,1),tt(i1-1:i1,2))
  if mod(i1,2)
    line(tt(i1,1),tt(i1,2),'Marker','o','LineStyle','none');
    delete(ht)
    % Note that sample number is shown as (i1-1)/2, since the state has
    % been saved after updating of a single parameter, instead of updating
    % all parameters before saving the state
    ht=text(-4,4,sprintf('Step %d',(i1-1)/2),'FontSize',16,'FontWeight','bold');
  end
  pause(2/(0.5+i1))
end

pause
% Illustrate warm-up
% first plot 
hb=line(tt(1,1),tt(1,2),'Marker','o','Color','m'); % this is to get nice legend
line(tt(1:2:50,1),tt(1:2:50,2),'Marker','o','LineStyle','none','Color','m');
line(tt(1:50,1),tt(1:50,2),'Color','m');
hl=legend([he h1 hb],'90% HPD','Markov chain','warm-up','location','southeast');
set(hl,'FontSize',14)
delete(ht)

pause
% Clear again
clf
% Plot 90% HPD
q=sort(sqrt(eig(S))*2.147);
he=ellipse(q(2),q(1),pi/4,0,0,'r'); 
% set axes
set(gca,'DataAspectRatio',[1 1 1])
set(gca,'XLim',[-4.5 4.5],'YLim',[-4.5 4.5])
set(gca,'Box','on')
% Plot again the first 200 samples withou warm-up
h1=line(tt(51:2:201,1),tt(51:2:201,2),'Marker','o','LineStyle','none');
title('Gibbs sampling')
hl=legend([he h1],'90% HPD','Samples from the chain after warm-up','location','southeast');
set(hl,'FontSize',14)

pause
% Clear again
clf
% Plot 90% HPD
q=sort(sqrt(eig(S))*2.147);
he=ellipse(q(2),q(1),pi/4,0,0,'r'); 
% set axes
set(gca,'DataAspectRatio',[1 1 1])
set(gca,'XLim',[-4.5 4.5],'YLim',[-4.5 4.5])
set(gca,'Box','on')
% Remove warm-up. In this case 50 first iterations, and plot the
% rest of the samples
h1=line(tt(51:2:end,1),tt(51:2:end,2),'Marker','o','LineStyle','none');
title('Gibbs sampling')
hl=legend([he h1],'90% HPD','950 samples from the chain','location','southeast');
set(hl,'FontSize',14)

pause
% Clear again
clf

subplot(3,1,1)
% show trends
burnin=51;
h1=plot([burnin:2:M]./2,tt(burnin:2:M,:));
line([burnin M]./2,[0 0])
axis tight
title('Gibbs sampling - trends')
hl=legend(h1,'\theta_1', '\theta_2');
set(hl,'FontSize',14)

subplot(3,1,2)
% show how the average behaves when more samples are obtained
burnin=51;
h1=plot([burnin:2:M]./2,cumsum(tt(burnin:2:M,:))./repmat([1:length(burnin:2:M)]',1,2));
line([burnin M]./2,[0 0])
axis tight
title('Cumulative average')
hl=legend(h1,'\theta_1', '\theta_2');
set(hl,'FontSize',14)

subplot(3,1,3)
% plot the estimated autocorrelation function
maxlag=20;
ac=acorr(tt(burnin:2:M,:),maxlag);
h1=plot(0:maxlag,[1 1;ac]);
set(gca,'XLim',[0 maxlag],'YLim',[-0.2 1])
line([0 maxlag],[0 0],'LineStyle','--','Color','k')
title('Estimate of the autocorrelation function')
hl=legend(h1,'\theta_1', '\theta_2');
set(hl,'FontSize',14)

pause
% Clear again
clf
% show how the average behaves when more samples are obtained
burnin=51;
h=plot([burnin:2:M]./2,cumsum(tt(burnin:2:M,:))./repmat([1:length(burnin:2:M)]',1,2));
line([burnin M]./2,[0 0])
axis tight
title('Gibbs sampling - cumulative average')
% compute the effective number of samples
[R,neff] = psrf(tt(burnin:2:end,:));
% for a simpler plot, use the minimum neff (each parameter has its own neff)
neff=min(neff);
% plot approximative 95% quantiles for Markov chain Monte Carlo error
% this takes into account that samples are not independent
line([[burnin:2:M]./2],1.96./sqrt([burnin:2:M]/2/(M./2./neff)),'LineStyle','--','Color','k')
h2=line([[burnin:2:M]./2],-1.96./sqrt([burnin:2:M]/2/(M./2./neff)),'LineStyle','--','Color','k');
legend([h' h2],'\theta_1','\theta_2','95% interval for MCMC error')
pause
% For comparison plot 95% quantiles for Monte Carlo error
% if we would be able to obtain independent samples
line([[burnin:2:M]./2],1.96./sqrt([burnin:2:M]/2),'LineStyle','-.','Color','k')
h3=line([[burnin:2:M]./2],-1.96./sqrt([burnin:2:M]/2),'LineStyle','-.','Color','k');
legend([h' h2 h3],'\theta_1','\theta_2','95% interval for MCMC error','95% interval for independent MC')
