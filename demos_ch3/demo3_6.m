% Bayesian data analysis
% Aki Vehtari <Aki.Vehtari@aalto.fi>

% Illustrate the posterior inference for Bioassay data (BDA3 p. 74-)

% data
x=[-.86 -.30 -.05 .73]';
n=[5 5 5 5]';
y=[0 1 3 5]';

clf
% default settings for the figure
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'Position',[620 60 400 640]);

% compute the posterior density in grid
%  - usually should be computed in logarithms!
%  - with alternative prior, check that range and spacing of A and B
%    are sensible
A=linspace(-4,8,50);
B=linspace(-10,40,50);
% for all combinations of A_i, B_j, i=1,...,200, j=1,...,200
for i=1:50
  a=A(i);
  for j=1:50
    b=B(j);
    p(j,i)=prod((1./(exp(-(a+b.*x))+1)).^y.*(1-1./(exp(-(a+b.*x))+1)).^(n-y));
  end
end

% plot the distribution in color
subplot(3,1,1)
pcolor(A,B,p);
shading interp
set(gcf,'renderer','painters')
xlim([-2 8])
ylim([-2 40])
xlabel('alpha')
ylabel('beta')

% plot samples from the distribution
subplot(3,1,2)
% catrand performs 2D-grid-sampling given values p in the grid
r=catrand(p,1000,1);
% r has indeces to a matrix (grid) length(A) x length(B)
% transform these to row and column indeces which can be used
% to index vectors A and B
[I,J]=ind2sub(size(p),r);
a=A(J);b=B(I);
% compute the grid spacing
ad=A(2)-A(1);bd=B(2)-B(1);
% add random jitter (see BDA3 p. 76)
plot(a+rand(size(a))*ad-ad/2,b+rand(size(b))*bd-bd/2,'.')
xlim([-2 8])
ylim([-2 40])
xlabel('alpha')
ylabel('beta')

pause
% plot the histogram of LD50 conditonal beta>0
subplot(3,1,3)
% logical vector for beta>0
bpi=b>0;
hist(-a(bpi)./b(bpi),[-0.5:0.02:0.5])
set(gca,'YTick',[],'XTick',[-0.4:0.2:0.4],'XLim',[-0.5 0.5])
xlabel('LD50 = -alpha/beta')

% Instructions for the exercise 3.11 in BDA3
% - Check that the range and spacing of A and B are sensible for the 
%   alternative prior
% - Compute the log-posterior in a grid
% - Scale the log-posterior by subtracting its maximum value before
%   exponentiating (see BDA3 p. 261)
% - Exponentiate
% - Normalize the posterior
% - Use 2D grid sampling 
% - In addition to the plots, report p(beta>0|x,y)
