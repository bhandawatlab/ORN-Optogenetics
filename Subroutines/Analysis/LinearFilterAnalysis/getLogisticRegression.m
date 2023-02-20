function [b,yhat] = getLogisticRegression(x,y,y_std,delay,fs)
% getLogisticRegression  calculates a logistic GLM fit. Note that y_std and
%   delay are only used for plotting purposes
%
%   Inputs: x = regressors (e.g. firing rate)
%           y = response (e.g. kinematics)
%           y_std = standard dev of output
%           delay = time period
%
%   Output: b = fitten parameters
%           yhat = predicted response
%   

[b,dev,stats] = glmfit(x,y,'binomial','link','logit');
yhat = glmval(b,x,'logit');

plot((-delay:delay)./fs,y,'-k','LineWidth',2);hold on;
plot((-delay:delay)./fs,yhat,'-g','LineWidth',2);
plot((-delay:delay)./fs,y-y_std,'--k','LineWidth',1);
plot((-delay:delay)./fs,y+y_std,'--k','LineWidth',1);
hold off;
legend({'empirical','logistic GLM fit'})
end