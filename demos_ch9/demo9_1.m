%   Author: Aki Vehtari <Aki.Vehtari@aalto.fi>

% Prof Gelman has a jar of coins. He promises that if the students
% guess how many coins there are, they will get all the coins in
% the jar. Students discuss and guess different values. Based on these
% they eventually present their uncertainty about the number of coins as
% a normal distribution N(160,30). What value they should guess?
x=50:260;
s=40;
px=normpdf(x,160,s);

clf
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold')

subplot(3,1,1)
% Plot the distribution
line([x;x],[zeros(size(x)); px],'color','b');
set(gca,'Ytick',[])
xlabel('Number of coins')
ylabel('Probability')
xlim([70 250])
ylim([0 0.018])

pause
% If students just want to guess right, and they do not care how
% much money tehy'll get they should guess the most probable value
h2=line([160 160],[0 normpdf(160,160,s)],'color','r');
legend(h2,'Most probable value = 160')

pause
% Alternatively students might want to maximize the exepected utility
% of the number coins. Assume that utility of the money is linear
subplot(3,1,2)
% Plot the utility 
line([x;x],[zeros(size(x)); x],'color','b');
xlim([70 250])
ylim([0 250])
xlabel('Guess')
ylabel('Utility if guess is correct')

pause
% If students guess value a, given their estimate of the
% uncertainity, probability that they get a coins is p(a),
% and expected utility is a*p(a)
subplot(3,1,3)
% Plot the expected utility 
line([x;x],[zeros(size(x)); x.*px],'color','b');
xlim([70 250])
ylim([0 3])
xlabel('Guess')
ylabel('Expected utility')

pause
% Compute the maximum of the expected utility
[meu,mi]=max(x.*px);
meux=x(mi);
%
h3=line([meux meux],[0 meu],'color','r');
legend(h3,sprintf('The guess which maximises the expected utility = %d',meux))
