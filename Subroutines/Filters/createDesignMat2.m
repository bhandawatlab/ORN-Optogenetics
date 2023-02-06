function [XStimAll,ntfilt] = createDesignMat2(stim,ntfilt,opts)
% createDesignMat2  Creates a design matrix X to use in solving for the 
%   linear filter b in the form y = Xb
%
%   Inputs: stim = n x time matrix of stimulus intensities
%           ntfilt = duration of filter
%           opts.ncells = number of flies (n)
%
%   Output: XStimAll = design matrix
%           ntfilt = duration of filter
%   

for cellnum = 1:opts.ncells
    Stim = stim(cellnum,:)';lenStim = length(Stim);
    lenStim2 = length(Stim)-ntfilt-1;
    
    paddedStim = [zeros(ntfilt,1); Stim; zeros(ntfilt,1)];
    c = 1:1:lenStim;
    r = lenStim:lenStim+2*ntfilt;
    tmp = hankel(c,r);
    Xstim = paddedStim(tmp);
    
    XStimAll((cellnum-1)*lenStim2+1:cellnum*lenStim2,:) = Xstim(1:end-ntfilt-1,:);
end
XStimAll = [ones(size(XStimAll,1),1),XStimAll];

end
