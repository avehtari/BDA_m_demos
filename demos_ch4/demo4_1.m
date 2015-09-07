% Bayesian data analysis
% Aki Vehtari <Aki.Vehtari@aalto.fi>

% Bioassay data, (BDA3 page 86)
x=[-.86 -.30 -.05 .73]';
n=[5 5 5 5]';
y=[0 1 3 5]';

clf
% default settings for the figure
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'Position',[200 60 820 640]);

% compute the posterior density in a grid
%  - usually this kind of computations should be computed in logarithms (see BDA3, p. 261)
%  - with alternative prior, check that range and spacing of A and B
%    are sensible
A=linspace(-4,8,200);
B=linspace(-10,40,200);
% for all combinations of A_i, B_j, i=1,...,200, j=1,...,200
for i=1:200
  a=A(i);
  for j=1:200
    b=B(j);
    p(j,i)=prod((1./(exp(-(a+b.*x))+1)).^y.*(1-1./(exp(-(a+b.*x))+1)).^(n-y));
  end
end

% plot the distribution in color
subplot(2,3,1)
pcolor(A,B,p);
shading interp
set(gcf,'renderer','painters')
axis([-2 6 -10 30])
xlabel('alpha')
ylabel('beta')

% plot samples from the distribution
subplot(2,3,2)
r=catrand(p,1000,1);
% r has indeces to a matrix (grid) length(A)xlength(B)
% transform these to row and column indeces which can be used
% index vectors A and B
[I,J]=ind2sub(size(p),r);
a=A(J);b=B(I);
% compute the grid spacing
ad=A(2)-A(1);bd=B(2)-B(1);
% add random jitter (see BDA3 76)
plot(a+rand(size(a))*ad-ad/2,b+rand(size(b))*bd-bd/2,'.')
axis([-2 6 -10 30])
xlabel('alpha')
ylabel('beta')
% beta describes how much increase in log dose has effect
% if beta were negative chemical would have healing properties!
% since none of the posterior sampes has negative beta,
% probability of negative beta is approximate to be less than 1/1000
% and correspondingly p(beta>0)>0.999
text(0,-7,'p(beta>0)>0.999','FontSize',14,'FontWeight','bold');

% plot histogram of LD50
subplot(2,3,3)
hist(-a./b,[-1:0.04:1])
set(gca,'YTick',[],'XTick',[-0.8:0.4:0.8],'XLim',[-0.8 0.8])
xlabel('LD50 = -alpha/beta')
drawnow


w0=[0 0];
% Compute gradients and Hessian analytically, and use Newton's
% method for optimisation
% You may use optimisation routines below for checking your results
% Uncomment desired lines
% For example, currently uncommented line provides, numerical approximation
% for the mode and Hessian, which you can compare to your results
% - no gradients provided (can't use LargeSacle)
opt=optimset('LargeScale','off');
% - if gradients are provided use following
%opt=optimset('GradObj','on');
% - if gradients and Hessian provided use following
%opt=optimset('GradObj','on','Hessian','on');
% - if gradients are provided and you want to check your gradients use following
%   (DerivativeCheck is not allowed with LargeScale)
%opt=optimset(opt,'LargeScale','off','DerivativeCheck','on');
% - optimize and get also Hessian H evaluated at optimum
[w,fval,exitflag,output,g,H]=fminunc(@bioassayfun,w0,opt,x,y,n);
% If using LargeScale without Hessian given, Hessian computed is sparse
if issparse(H)
  H=full(H);
end
S=inv(H);

% compute the normal approximation density in grid
% this is just for the illustration
A=linspace(-4,8,200);
B=linspace(-10,40,200);
% for all combinations of A_i, B_j, i=1,...,200, j=1,...,200
for i=1:200
  a=A(i);
  for j=1:200
    b=B(j);
    p(j,i)=exp(mnorm_lpdf([a b],w,S));
  end
end

% plot the distribution in color
subplot(2,3,4)
pcolor(A,B,p);
shading interp
set(gcf,'renderer','painters')
axis([-2 6 -10 30])
xlabel('alpha')
ylabel('beta')

% plot samples from the distribution
subplot(2,3,5)
r=mvnrnd(w,S,1000);
a=r(:,1);b=r(:,2);
plot(a,b,'.')
axis([-2 6 -10 30])
xlabel('alpha')
ylabel('beta')
% Normal approximation does not take into account that the posterior
% is not symmetric and that there is very low density for negative
% beta values. Based on the draws from the normal approximation
% is is estimated that there is about 6% probability that beta is negative!
text(0,-7,sprintf('p(beta>0)=%.2f',mean(b>0)),'FontSize',14,'FontWeight','bold');

% plot histogram of LD50
subplot(2,3,6)
% Since we have strong prior belief that beta should not be negative we can
% improve our normal approximation by conditioning on beta>0
% logical vector for beta>0
bpi=b>0;
hist(-a(bpi)./b(bpi),[-1:0.04:1])
set(gca,'YTick',[],'XTick',[-0.8:0.4:0.8],'XLim',[-0.8 0.8])
xlabel(['LD50 = -alpha/beta';'            beta>0'])

