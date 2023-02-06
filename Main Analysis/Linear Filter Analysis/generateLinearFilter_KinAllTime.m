function [b] = generateLinearFilter_KinAllTime(stim_all,kin,kin_std,ntfilt,filter_type,yl,optFitN,U,s,V,XStimNew,fs)
% generateLinearFilter_KinAllTime  calculates the linear filter from firing
%   rate to kinematics
%
%   Inputs: stim_all = input (i.e. firing rate)
%           kin = output (i.e. speed,curv,etc)
%           kin_std = standard deviation of output (plotting purposes)
%           ntfilt = filter duration (in frames)
%           filter_type = label for type of filter (plotting purposes)
%           yl = ylabel (plotting purposes)
%           optFitN = user defined lambda value (unused)
%           U,s,V = singular value decomposition matrices
%           XStimNew = filter design matrix
%           fs = sampling rate
%
%   Output: b = linear filter for each lambda value
%   

% solve for the ridge regression
nLamb = 30;
lamb = 2.^((1:nLamb)-15);
ntfilt_all = numel(s);
b = zeros(ntfilt_all,nLamb);rho = zeros(1,nLamb);eta = zeros(1,nLamb);
for i = 1:nLamb
    [b(:,i),rho(i),eta(i)] = tikhonov(U,s,V,kin,lamb(i),zeros(ntfilt_all,1));
end
% choose some lambdas to look at
optFitN_all = 2:6:20;
figure;set(gcf,'Position',[2 42 838 924]);
for fitN = 1:numel(optFitN_all)
    b2 = b(:,optFitN_all(fitN));
    % apply the filter to the data and reshape to nTrials x nPts matrix
    turnProb_fit = XStimNew*b2;
    
    b_broken = reshape(b2(2:end),[],numel(stim_all));
    
    subplot(5,2,(fitN-1).*2+1);
    plot([-ntfilt:1:-1]./fs,b_broken);
    xlabel('delay (s)');
    title(['lambda=' num2str(lamb(optFitN_all(fitN)))])
    subplot(5,2,fitN.*2);
    plot((1:numel(kin))./fs,kin,'-k','LineWidth',2);hold on;
    plot((1:numel(kin))./fs,turnProb_fit,'-g','LineWidth',2);
    plot((1:numel(kin))./fs,kin+kin_std,'--k','LineWidth',1);
    plot((1:numel(kin))./fs,kin-kin_std,'--k','LineWidth',1);
    ylabel(yl);
    xlabel('delay (s)');
    title(['err = ' num2str(rho(optFitN_all(fitN)))])
    
    if fitN == 1
        subplot(5,2,(fitN-1).*2+1);
        legend(filter_type)
        subplot(5,2,fitN.*2);
        legend({'emp','fit'})
    end
end
if isempty(optFitN)
    optFitN = 16;
end
lambLabel = cellfun(@(x) strtrim(x), cellstr(num2str(2.^(optFitN_all-15)')), 'UniformOutput', false);

subplot(5,2,9)
loglog(rho,eta);hold on;
scatter(rho(optFitN_all),eta(optFitN_all),'ro');
text(rho(optFitN_all),eta(optFitN_all)./5,lambLabel,'HorizontalAlignment','center');
xlabel('residual norm || A x - b ||_2')
ylabel('solution norm || x ||_2')
title([sprintf('fitting with: %s,  %s ',filter_type{:})])

end