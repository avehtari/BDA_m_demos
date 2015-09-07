% Bayesian data analysis
% Aki Vehtari <Aki.Vehtari@aalto.fi>

% Hierarchical model for Rats experiment (BDA3, p. 102)

% Rat data (BDA3, p. 102)
y=[0  0  0  0  0  0  0  0  0  0 ...
0  0  0  0  1  1  1  1  1  1 ...
1  1  2  2  2  2  2  2  2  2 ...
2  1  5  2  5  3  2  7  7  3 ...
3  2  9  10  4  4  4  4  4  4 ...
4  10  4  4  4  5  11  12  5  5 ...
6  5  6  6  6  6  16  15  15  9 ...
4];
n=[20  20  20  20  20  20  20  19  19  19 ...
19  18  18  17  20  20  20  20  19  19 ...
18  18  25  24  23  20  20  20  20  20 ...
20  10  49  19  46  27  17  49  47  20 ...
20  13  48  50  20  20  20  20  20  20 ...
20  48  19  19  19  22  46  49  20  20 ...
23  19  22  20  20  20  52  46  47  24 ...
14];

% Default setting for figure
clf
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'Position',[200 60 820 640]);
box on

% How many experiments
M=length(y);

clf
x=linspace(0,1,1000);
% Plot separate model
box on
for i=1:M
  px=betapdf(x,y(i)+1,n(i)-y(i)+1);
  h1=line(x,px);
end
% The last one is for emphasize colored red
h2=line(x,px,'Color','r');
set(gca,'Ytick',[])
title('Tumor rats - Separate model')
legend([h1 h2],'Posterior of \theta_j','Posterior of \theta_{71}')
xlabel('\theta')

pause
% Plot both separate and pooled model
% First separate model (same as above)
subplot(2,1,1)
box on
for i=1:M
  px=betapdf(x,y(i)+1,n(i)-y(i)+1);
  h1=line(x,px);
end
% The last one is for emphasize colored red
h2=line(x,px,'Color','r');
set(gca,'Ytick',[])
title('Separate model')
legend([h1 h2],'Posterior of \theta_j','Posterior of \theta_{71}')
xlabel('\theta')
% Then pooled model
subplot(2,1,2)
px=betapdf(x,sum(y)+1,sum(n)-sum(y)+1);
% Common theta for all
line(x,px,'Color','r');
set(gca,'Ytick',[])
title('Pooled model')
legend('Posterior of common \theta')
xlabel('\theta')
drawnow

% Compute the marginal posterior of alpha and beta in hierarchical model
% Use grid 
A=linspace(0.5,6,100);
B=linspace(3,33,100);
for i=1:100
  a=A(i);
  for j=1:100
    b=B(j);
    % Use logarithms for numerical accuracy!
    lp(j,i)=log(a+b).*(-5/2)+sum(gammaln(a+b)-gammaln(a)-gammaln(b)+gammaln(a+y)+gammaln(b+n-y)-gammaln(a+b+n));
  end
end
% Subtract maximum value to avoid over/underflow in exponentation
lp=lp-max(lp(:));
p=exp(lp);

% Notify that previous computation has been completed
fprintf('.')
% and then wait for enter
pause
% Plot the marginal posterior
clf
pcolor(A,B,p);
shading interp
set(gcf,'renderer','painters')
set(gca,'Xtick',1:6)
xlabel('alpha')
ylabel('beta')
title('The marginal posterior of alpha and beta in hierarchical model')

% Plot samples from the distribution
% catrand performs 2D-grid-sampling given values p in the grid
r=catrand(p,1000,1);
% r has indeces to a matrix (grid) length(A)xlength(B)
% transform these to row and column indeces which can be used
% index vectors A and B
[I,J]=ind2sub(size(p),r);
a=A(J);b=B(I);
pause
clf
% Plot samples from the distribution of distributions Beta(alpha,beta),
% that is, plot Beta(alpha,beta) using posterior samples of alpha and beta
subplot(2,1,1)
box on
x=linspace(0,1,1000);
for i=1:20
  px=betapdf(x,a(i),b(i));
  line(x,px)
end
set(gca,'Ytick',[])
title('Posterior samples from the distribution of distributions Beta(\alpha,\beta)')

% The average of above distributions, is the predictive distribution
% for a new theta, and also the prior distribution for theta_j
% Plot this
subplot(2,1,2)
box on
x=linspace(0,1,1000);
pxs=zeros(100,1000);
for i=1:100
  pxs(i,:)=betapdf(x,a(i),b(i));
end
px=mean(pxs);
line(x,px)
set(gca,'Ytick',[])
title(['Predictive distribution for a new \theta and prior for \theta_j'])

pause

% And finally compare the separate model and hierarchical model
% First plot the separate model (same as above)
clf
subplot(2,1,1)
box on
x=linspace(0,1,1000);
pxs=zeros(100,1000);
% Note that for clarity only every 7th distribution is plotted
for i=8:7:M
  px=betapdf(x,y(i)+1,n(i)-y(i)+1);
  line(x,px)
end
line(x,px,'Color','r')
xlim([0 1])
set(gca,'Ytick',[])
title('Separate model')

% And the hierarchical model 
% Note that these marginal posteriors for theta_j are more narrow than
% in the separate model case, due to the borrowed information from the
% other theta_j's
subplot(2,1,2)
box on
x=linspace(0,0.6,200);
pxs=zeros(100,200);
% Note that for clarity only every 7th distribution is plotted
for i=8:7:M
  for j=1:100
    pxs(j,:)=betapdf(x,y(i)+a(j),n(i)-y(i)+b(j));
  end
  px=mean(pxs);
  line(x,px)
end
line(x,px,'Color','r')
set(gca,'Ytick',[])
xlim([0 1])
title('Hierarchical model')
