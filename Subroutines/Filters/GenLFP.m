function [LFP_fit,ntfilt] = GenLFP(stim,b,ntfilt,ncells,fs)
% get the LFP
opts.ncells = ncells;
if isempty(ntfilt)
    ntfilt = 5.*fs;
end
nthist = ntfilt;%nthist;
sps = zeros(size(stim'));

% choose the first 1.5 s for the filter
b2 = b(:,end);%b([1 10:end]);%b([1,(end-emp.fs1*6+1):end]);
[XStimAll,~,~,~] = createDesignMat(stim,sps,ntfilt,nthist,opts);
XStimNew = XStimAll;%XStimAll(:,[1 10:end]);%XStimAll(:,(end-emp.fs1*6+1):end);

% apply the filter to the data
% yfit = glmval(b2,XStimNew,'identity');
yfit = XStimNew*b2;
yfit = reshape(yfit,[],ncells);
LFP_fit = yfit;

end








