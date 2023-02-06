function [XStimAll,XdsgnAll,ntfilt,nthist] = createDesignMat(stim,resp,ntfilt,nthist,opts)
% createDesignMat2  Creates a design matrix X to use in solving for the 
%   linear filter b in the form y = Xb. Includes the option of history
%
%   Inputs: stim = n x time matrix of stimulus intensities
%           resp = n x time matrix of response
%           ntfilt = duration of filter
%           nthist = duration of response filter
%           opts.ncells = number of flies (n)
%
%   Output: XStimAll = design matrix without response history
%           XdsgnAll = design matrix with response history
%           ntfilt = duration of filter
%           nthist = duration of response filter
%   

for cellnum = 1:opts.ncells
    Stim = stim(cellnum,:)';lenStim = length(Stim);
    paddedStim = [ones(ntfilt-1,1).*Stim(1); Stim]; % pad early bins of stimulus with the first value
    Xstim = hankel(paddedStim(1:end-ntfilt+1), Stim(end-ntfilt+1:end));
    
    XStimAll((cellnum-1)*lenStim+1:cellnum*lenStim,:) = Xstim;
    
    if ~isempty(resp)
        % Build spike-history design matrix
        paddedSps = [zeros(nthist,1); resp(1:end-1,cellnum)];
        Xsp = hankel(paddedSps(1:end-nthist+1), paddedSps(end-nthist+1:end));
        
        % Combine these into a single design matrix
        Xdsgn = [Xstim,Xsp];
        XdsgnAll((cellnum-1)*lenStim+1:cellnum*lenStim,:) = Xdsgn;
        
    end
    
end
XStimAll = [ones(size(XStimAll,1),1),XStimAll];
if ~isempty(resp)
    XdsgnAll = [ones(size(XdsgnAll,1),1),XdsgnAll];
else
    XdsgnAll = [];
end

end