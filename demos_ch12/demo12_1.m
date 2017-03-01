% Bayesian data analysis
% Aki Vehtari <Aki.Vehtari@aalto.fi>

% HMC/NUTS demonstration

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
tt=zeros(M,2);   % variable for saving samples
tt=[t1 t2];      % save starting point

% Hamiltonian Monte Carlo %
% Using hmc2 function from GPstuff 
% http://research.cs.aalto.fi/pml/software/gpstuff/
func=@(x) -mnorm_lpdf(x,[y1 y2],S); % negative log density
grad=@(x) (S\x')';                  % gradient
% hmc options
% note that we are calling hmc2 with unusual options
% usually steps>1 and decay<1, but because the intermediate steps are
% not stored by hmc2-function, here steps=1 and decay=1 to be able to
% visualise each step in HMC
opt=hmc2_opt;
opt.steps=1;     
opt.stepadj=0.05;
opt.nsamples=40;
opt.persistence=1;
opt.decay=1;
fprintf('Running HMC')
for i1=1:1000
    [samples, energies, diagn]=hmc2(func, tt(end,:), opt, grad);
    tt=[tt;samples];
    % usually energy would be sampled after a number of steps when
    % persistence=0 or decay<0, but due to special options for
    % visualisation the moment is manually reseted here
    hmc_state=hmc2('state');
    hmc_state.mom=[];
    hmc2('state',hmc_state);
end
fprintf(' - done.\n')

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
i1=1;
hs=line(tt(1,1),tt(1,2),'Marker','o','LineStyle','none');
title('Hamiltonian Monte Carlo')
xlabel('\theta_1')
ylabel('\theta_2')
hl=legend([he hs],'90% HPD','Starting point','location','southeast');
set(hl,'FontSize',14)
% wait for enter
pause

% plot progression 1:40
h1=line(tt(i1,1),tt(i1,2),'LineStyle','-','Marker','none','color',[0 0.5 0]);
h2=line(tt(i1,1),tt(i1,2),'LineStyle','none','Marker','o');
hl=legend([he hs],'90% HPD','Samples','location','southeast');
for i1=1:40
    % plot the trajectory
    set(h1,'XData',tt(max(1,i1-100):i1,1),'YData',tt(max(1,i1-100):i1,2))
    if i1>=40
        set(h2,'XData',tt(40:40:i1,1),'YData',tt(40:40:i1,2))
    end
    pause(0.1)
end
% wait for enter
pause

% plot progression 41:80
for i1=41:80
    % plot the trajectory
    set(h1,'XData',tt(max(1,i1-100):i1,1),'YData',tt(max(1,i1-100):i1,2))
    if i1>=40
        set(h2,'XData',tt(40:40:i1,1),'YData',tt(40:40:i1,2))
    end
    pause(0.1)
end
% wait for enter
pause

% plot progression 81:120
for i1=81:120
    % plot the trajectory
    set(h1,'XData',tt(max(1,i1-100):i1,1),'YData',tt(max(1,i1-100):i1,2))
    if i1>=40
        set(h2,'XData',tt(40:40:i1,1),'YData',tt(40:40:i1,2))
    end
    pause(0.1)
end
% wait for enter
pause

% plot progression 121:1000
for i1=121:1000
    % plot the trajectory
    set(h1,'XData',tt(max(1,i1-100):i1,1),'YData',tt(max(1,i1-100):i1,2))
    if i1>=40
        set(h2,'XData',tt(40:40:i1,1),'YData',tt(40:40:i1,2))
    end
    pause(0.01)
end
% wait for enter
pause
% delete trajectory lines
delete(h1)

% wait for enter
pause

% plot all iterations without lines
delete(hs)
set(h2,'XData',tt(100:20:end,1),'YData',tt(100:20:end,2))
hl=legend([he h2],'90% HPD','Samples after warm-up','location','southeast');
% wait for enter
pause

% No-U-Turn-Sampling
tt=[t1 t2]; % start again from the initial starting point
fg=@(x) deal(mnorm_lpdf(x,[y1 y2],S), -(S\x')'); % log density and gradient
% options
% Using hmc_nuts function from GPstuff 
% http://research.cs.aalto.fi/pml/software/gpstuff/
% the algorithm parameters are automatically adapted
% and hmc_nuts is called in a usual way as we don't visualise the
% intermediate steps
opt.M=M;
opt.Madapt=100;
fprintf('\nRunning NUTS')
[samples, energies, diagn]=hmc_nuts(fg, tt(1,:), opt);
fprintf(' - done.\n')
tt=samples;

% plot NUTS draws
subplot(3,1,1)
% show trends
burnin=105;
h1=plot([burnin:M],tt(burnin:M,:));
line([burnin M],[0 0])
axis tight
title('HMC-NUTS - trends')
hl=legend(h1,'\theta_1', '\theta_2');
set(hl,'FontSize',14)

subplot(3,1,2)
% show how the average behaves when more samples are obtained
burnin=105;
h1=plot([burnin:M],cumsum(tt(burnin:M,:))./repmat([1:length(burnin:M)]',1,2));
line([burnin M],[0 0])
axis([105 5000 -0.5 0.2])
title('Cumulative average')
hl=legend(h1,'\theta_1', '\theta_2','location','southeast');
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

% wait for enter
pause
clf
% show how the average behaves when more samples are obtained
burnin=105;
h=plot([burnin:M],cumsum(tt(burnin:M,:))./repmat([1:length(burnin:M)]',1,2));
line([burnin M],[0 0])
axis([105 5000 -0.5 0.5])
legend('\theta_1','\theta_2')
title('HMC-NUTS sampling - cumulative average')
% compute the effective number of samples
[R,neff] = psrf(tt(burnin:M,:));
% for a simpler plot, use the minimum neff (each parameter has its own neff)
neff=min(neff);
% plot approximative 95% quantiles for Markov chain Monte Carlo error
% this takes into account that samples are not independent
line([burnin:M],1.96./sqrt([burnin:M]/(M/neff)),'LineStyle','--','Color','k')
h2=line([burnin:M],-1.96./sqrt([burnin:M]/(M/neff)),'LineStyle','--','Color','k');
legend([h' h2],'\theta_1','\theta_2','95% interval for MCMC error')
% wait for enter
pause
% For comparison plot 95% quantiles for Monte Carlo error
% if we would be able to obtain independent samples
line([burnin:M],1.96./sqrt([burnin:M]),'LineStyle','-.','Color','k')
h3=line([burnin:M],-1.96./sqrt([burnin:M]),'LineStyle','-.','Color','k');
legend([h' h2 h3],'\theta_1','\theta_2','95% interval for MCMC error','95% interval for independent MC')
