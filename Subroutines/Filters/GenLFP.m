function [LFP_fit,ntfilt] = GenLFP(stim,b,ntfilt,nflys,fs)
% GenLFP  Computes the local field potential (LFP) from a stimulus and a
%   LFP filter
%
%   Inputs: stim = n x time matrix of stimulus intensities
%           b = tau+1 x 1 vector of linear filter weights (DC offset = b(1))
%           ntfilt = duration of filter (tau)
%           ncells = number of flies (n)
%           fs = sampling rate
%
%   Output: LFP_fit = time x n matrix of predicted LFPs
%           ntfilt = duration of filter (tau)
%   

% get the LFP
opts.ncells = nflys;
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
yfit = reshape(yfit,[],nflys);
LFP_fit = yfit;

end








