function [XStimAll,XdsgnAll,ntfilt,nthist] = createDesignMat(stim,sps,ntfilt,nthist,opts)

for cellnum = 1:opts.ncells
    Stim = stim(cellnum,:)';lenStim = length(Stim);
    % Build stimulus design matrix (using 'hankel');
    paddedStim = [ones(ntfilt-1,1).*Stim(1); Stim]; % pad early bins of stimulus with the first value
    Xstim = hankel(paddedStim(1:end-ntfilt+1), Stim(end-ntfilt+1:end));
    
    XStimAll((cellnum-1)*lenStim+1:cellnum*lenStim,:) = Xstim;
    
    if ~isempty(sps)
        % Build spike-history design matrix
        paddedSps = [zeros(nthist,1); sps(1:end-1,cellnum)];
        % SUPER important: note that this doesn't include the spike count for the
        % bin we're predicting? The spike train is shifted by one bin (back in
        % time) relative to the stimulus design matrix
        Xsp = hankel(paddedSps(1:end-nthist+1), paddedSps(end-nthist+1:end));
        
        % Combine these into a single design matrix
        Xdsgn = [Xstim,Xsp];
        XdsgnAll((cellnum-1)*lenStim+1:cellnum*lenStim,:) = Xdsgn;
        
    end
    
end
XStimAll = [ones(size(XStimAll,1),1),XStimAll];
if ~isempty(sps)
    XdsgnAll = [ones(size(XdsgnAll,1),1),XdsgnAll];
else
    XdsgnAll = [];
end

end