function [sps_pred] = genFiringRate(LFP_fit,b2cons,ncells,ntFilt)
% genFiringRate  Computes the firing rate from a LFP and a firing rate
%	filter
%
%   Inputs: LFP_fit = time x n matrix of LFPs
%           b2cons = tau+1 x 1 vector of linear filter weights 
%               (DC offset = b(1))
%           ntfilt = duration of filter (tau)
%           ncells = number of flies (n)
%
%   Output: sps_pred = n x time matrix of predicted firing rates
%   

opts.ncells = ncells;

% design matrix
[XStimAll,~,~,~] = createDesignMat(LFP_fit',[],ntFilt,ntFilt,opts);

sps_pred = XStimAll*b2cons;
sps_pred = reshape(sps_pred,[],ncells)';
sps_pred(sps_pred<0) = 0;

end