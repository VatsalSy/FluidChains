%% This coode reads the eta variable for all the cases and developes a corelation
%% model
clc
clear
close
load('eta')
params = sortrows(params);
%remove = params(213:218,3:3);
Xp = X;
%remove = [remove;22];
Xp(22,:) = [];
etap = eta;
etap(22,:) = [];

etan = log(etap(:,1));

Xt = Xp;
Xt(:,1) = log(sin(pi/180*Xt(:,1)));
Xt(:,2) = log(Xt(:,2));
Xt(:,3) = log(Xt(:,3));
Xt(:,4) = log(Xt(:,4));
temp = [ones(size(Xt,1),1) Xt];
Xt = temp;
[b1,~,~,~,stats1] = regress(etan,Xt);
save('finalRegress.mat','b1','X','eta'); % saves the regression analysis results of this code in finalRegress.mat
stats1
etaTh = zeros(size(X,1),1);
for i = 1:1:size(X,1)
    etaTh(i) = exp(b1(1))*sin(pi/180*X(i,1))^b1(2)*X(i,2)^b1(3)*X(i,3)^b1(4)*X(i,4)^b1(5);
end
plot(etaTh,eta(:,1),'k.','MarkerSize',10)
hold on
x = linspace(min(etaTh),max(etaTh),1200);
plot(x,x,'k-','LineWidth',3)
plot(x,1.1*x,'k-','LineWidth',3)
plot(x,0.9*x,'k-','LineWidth',3)
% %b1

% subplot(2,2,1)
% plot(Xt(:,2),etan,'k*')
% subplot(2,2,2)
% plot(Xt(:,3),etan,'k*')
% subplot(2,2,3)
% plot(Xt(:,3),etan,'k*')
% subplot(2,2,4)
% plot(Xt(:,4),etan,'k*')
