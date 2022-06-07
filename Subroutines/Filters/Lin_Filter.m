function [sps_pred] = Lin_Filter(LFP_fit,R,b2cons,ncells,ntFilt,fs)

%Stim = LFP_fit';
Stim = [LFP_fit;repmat(LFP_fit(end,:),ntFilt+1,1)]';
opts.ncells = ncells;

% design matrix
%[XStimAll,~] = createDesignMat2(Stim,ntFilt,opts);
[XStimAll,~,~,~] = createDesignMat(LFP_fit',[],ntFilt,ntFilt,opts);

sps_pred = XStimAll*b2cons;
sps_pred = reshape(sps_pred,[],ncells)';
sps_pred(sps_pred<0) = 0;

% nT = size(R,2); % number of time bins in stimulus
% dtStim = 1./fs;
% tbins = (.5:nT)*dtStim;
% figure;plot([-ntFilt:1:size(XStimAll,2)-ntFilt-2].*dtStim,b2cons(2:end));title(['y-intercept=' num2str(b2cons(1))])
% grid on
% figure;
% for cellnum = 1:6
%     subplot(3,2,cellnum)
%     yyaxis left
%     plot(tbins(1:size(sps_pred,2)),sps_pred(cellnum,:),'-','Color',[0, 0.4470, 0.7410]);hold on
%     plot(tbins,R(cellnum,:),'-k');ylim([0 60])
%     yyaxis right
%     plot(tbins,Stim(cellnum,1:size(R,2)));ylim([0 4])
% end

end