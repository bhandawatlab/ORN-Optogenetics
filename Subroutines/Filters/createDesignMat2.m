function [XStimAll,ntfilt] = createDesignMat2(stim,ntfilt,opts)

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
