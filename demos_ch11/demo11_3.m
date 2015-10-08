% Bayesian data analysis
% Aki Vehtari <Aki.Vehtari@aalto.fi>

% Convergence demonstration

clf
% default settings for figure
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')

y1=0;y2=0;r=0.8; % parameters of the Normal-distribution, which is used as a
                 % toy target distribution for illustration
S=[1 r;r 1];     % covariance matrix of the toy distribution

t1=-2.5;t2=2.5;  % starting value for the chain
sp=0.8;          % scale for proposal distribution
M=5000;          % number of iterations
tt=zeros(M,2);   % variable for saving samples
tt(1,:)=[t1 t2]; % save starting point
rr=0;            % variable for saving number of rejections

Y1=linspace(-4,4,200);
Y2=linspace(-4,4,200);

% different starting points
t1s=[-2.5 -2.5 2.5 2.5];
t2s=[2.5 -2.5 2.5 -2.5];
sp=0.2;
%for i2=1:4
%  tt=zeros(M,2);   % variable for saving samples
%  tt(1,:)=[t1s(i2) t2s(i2)]; % save starting point
%  rr=0;            % variable for saving number of rejections
%  % Metropolis algorithm here
%  tts(:,:,i2)=tt;
%end
% precomputed values from 4 Metropolis chains
% note that for the illustration purposes the proposal distribution was
% intentionally selected to get slow convergence
load demo11_3
%>> whos -file demo11_3
%  Name      Size                           Bytes  Class
%  tt     5000x2                            80000  double array
%  tts    5000x2x4                         320000  double array

plot(squeeze(tts(1:50,1,:)))
title('Multiple chains - no convergence yet')
xlabel('t')
pause

plot(squeeze(tts(1:1000,1,:)))
title('Multiple chains - probabaly converged')
xlabel('t')
pause

burnin=201;
h=plot(burnin:5000,cumsum(squeeze(tts(burnin:5000,1,:)))./repmat([1:(5000-burnin+1)]',1,4));
xlabel('t')
title('Multiple chains - estimates given t first samples after warm-up')

pause
line([burnin M],[0 0])
axis([0 5000 -1 1])
burnin=100;
% compute the effective number of samples
[R,neff] = psrf(tt(burnin:M,:));
% for a simpler plot, use the minimum neff (each parameter has its own neff)
neff=min(neff);
% plot approximative 95% quantiles for Markov chain Monte Carlo error
% this takes into account that samples are not independent
line([burnin:M],1.96./sqrt([burnin:M]/(M/neff)),'LineStyle','--','Color','k')
h2=line([burnin:M],-1.96./sqrt([burnin:M]/(M/neff)),'LineStyle','--','Color','k');
legend([h' h2],'1','2','3','4','95% interval for MCMC error')
xlabel('t')
